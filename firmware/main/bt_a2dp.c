/**
 * bt_a2dp.c — Bluetooth A2DP sink implementation
 *
 * Architecture:
 *   - Bluedroid stack runs on one core (internal BT task).
 *   - The SBC decoder calls our A2DP data callback from the BT task.
 *   - Data is pushed into a FreeRTOS StreamBuffer, which wakes a
 *     dedicated DSP+I2S task on the application core.
 *   - This decouples Bluetooth timing from DSP+I2S timing.
 */

#include "bt_a2dp.h"
#include "eq_coeffs.h"

#include <string.h>
#include "freertos/FreeRTOS.h"
#include "freertos/stream_buffer.h"
#include "freertos/task.h"

#include "esp_log.h"
#include "esp_timer.h"
#include "nvs_flash.h"
#include "esp_bt.h"
#include "esp_bt_main.h"
#include "esp_bt_device.h"
#include "esp_gap_bt_api.h"
#include "esp_a2dp_api.h"
#include "esp_avrc_api.h"

static const char *TAG = "bt_a2dp";

/* ── StreamBuffer between BT task and DSP task ────────────────────── */

/* Hold ~100ms of stereo 16-bit audio at 48 kHz (worst case).
 * 100ms × 48000 × 4 bytes/frame ≈ 19200 → round up to 20480.     */
#define PCM_STREAM_BUF_SIZE  20480
#define PCM_CHUNK            512    /* bytes per read (128 stereo frames) */

static StreamBufferHandle_t pcm_stream;
static bt_a2dp_data_cb_t    user_data_cb;
static bt_event_cb_t        user_event_cb;
static int                  current_sample_rate = 44100;
static bool                 audio_running = false;

/* ── AVRCP absolute volume (phone → sink) ─────────────────────── */

static volatile uint8_t s_volume = 127;  /* 0=mute, 127=max; phone overwrites via AVRCP */

/* ── Forward declarations ──────────────────────────────────────────── */

static void bt_avrc_tg_cb(esp_avrc_tg_cb_event_t event,
                          esp_avrc_tg_cb_param_t *param);
static void bt_a2dp_event_cb(esp_a2d_cb_event_t event,
                             esp_a2d_cb_param_t *param);
static void bt_a2dp_data_cb(const uint8_t *data, uint32_t len);
static void dsp_i2s_task(void *arg);

/* ───────────────────────────────────────────────────────────────────
 *  Public API
 * ─────────────────────────────────────────────────────────────────── */

void bt_a2dp_sink_init(const char *device_name, bt_a2dp_data_cb_t data_cb)
{
    user_data_cb = data_cb;

    /* ── NVS (needed by Bluetooth stack) ──────────────────────── */
    esp_err_t err = nvs_flash_init();
    if (err == ESP_ERR_NVS_NO_FREE_PAGES ||
        err == ESP_ERR_NVS_NEW_VERSION_FOUND) {
        nvs_flash_erase();
        nvs_flash_init();
    }

    /* ── Bluetooth controller + Bluedroid ─────────────────────── */
    esp_bt_controller_config_t bt_cfg = BT_CONTROLLER_INIT_CONFIG_DEFAULT();
    ESP_ERROR_CHECK(esp_bt_controller_init(&bt_cfg));
    ESP_ERROR_CHECK(esp_bt_controller_enable(ESP_BT_MODE_CLASSIC_BT));
    ESP_ERROR_CHECK(esp_bluedroid_init());
    ESP_ERROR_CHECK(esp_bluedroid_enable());

    /* ── Device name ──────────────────────────────────────────── */
    esp_bt_dev_set_device_name(device_name);

    /* ── AVRCP target (volume sync) ──────────────────────────── */
    /* Must be initialized before A2DP so the SDP record
     * advertises absolute volume support during pairing. */
    ESP_ERROR_CHECK(esp_avrc_tg_init());
    ESP_ERROR_CHECK(esp_avrc_tg_register_callback(bt_avrc_tg_cb));

    /* Advertise VOLUME_CHANGE notification capability.
     * Without this the phone doesn't know the sink supports
     * absolute volume → falls back to stream-relative attenuation
     * → enhancer starved of signal → fundamental leakage. */
    esp_avrc_rn_evt_cap_mask_t evt_cap = {0};
    esp_avrc_rn_evt_bit_mask_operation(ESP_AVRC_BIT_MASK_OP_SET, &evt_cap,
                                        ESP_AVRC_RN_VOLUME_CHANGE);
    ESP_ERROR_CHECK(esp_avrc_tg_set_rn_evt_cap(&evt_cap));

    /* ── A2DP sink ────────────────────────────────────────────── */
    ESP_ERROR_CHECK(esp_a2d_register_callback(bt_a2dp_event_cb));
    ESP_ERROR_CHECK(esp_a2d_sink_register_data_callback(bt_a2dp_data_cb));
    ESP_ERROR_CHECK(esp_a2d_sink_init());

    /* ── StreamBuffer + DSP task (runs on app core) ───────────── */
    pcm_stream = xStreamBufferCreate(PCM_STREAM_BUF_SIZE, PCM_CHUNK);
    xTaskCreatePinnedToCore(dsp_i2s_task, "dsp_i2s", 4096, NULL, 5,
                            NULL, 1);

    /* ── Set discoverable ─────────────────────────────────────── */
    esp_bt_gap_set_scan_mode(ESP_BT_CONNECTABLE,
                             ESP_BT_GENERAL_DISCOVERABLE);

    ESP_LOGI(TAG, "Ready as \"%s\". Waiting for connection…", device_name);
}

void bt_a2dp_set_event_callback(bt_event_cb_t cb)
{
    user_event_cb = cb;
}

uint8_t bt_a2dp_get_volume(void)
{
    return s_volume;
}

/* ───────────────────────────────────────────────────────────────────
 *  AVRCP target callback — runs in Bluedroid task context
 * ─────────────────────────────────────────────────────────────────── */

static void bt_avrc_tg_cb(esp_avrc_tg_cb_event_t event,
                          esp_avrc_tg_cb_param_t *param)
{
    switch (event) {
    case ESP_AVRC_TG_SET_ABSOLUTE_VOLUME_CMD_EVT:
        s_volume = param->set_abs_vol.volume;
        ESP_LOGI(TAG, "Volume set to %u/127 (%.0f%%)",
                 (unsigned)s_volume,
                 (double)s_volume * 100.0 / 127.0);
        break;

    case ESP_AVRC_TG_REGISTER_NOTIFICATION_EVT:
        /* Phone asks to be notified when our volume changes.
         * Respond immediately with our current volume so the
         * phone's slider matches on connect.  INTERIM means
         * "here's the current value, I'll tell you if it changes." */
        if (param->reg_ntf.event_id == ESP_AVRC_RN_VOLUME_CHANGE) {
            esp_avrc_rn_param_t rn = { .volume = s_volume };
            esp_avrc_tg_send_rn_rsp(ESP_AVRC_RN_VOLUME_CHANGE,
                                    ESP_AVRC_RN_RSP_INTERIM, &rn);
        }
        break;

    default:
        break;
    }
}

/* ───────────────────────────────────────────────────────────────────
 *  A2DP event handler
 * ─────────────────────────────────────────────────────────────────── */

static void bt_a2dp_event_cb(esp_a2d_cb_event_t event,
                             esp_a2d_cb_param_t *param)
{
    switch (event) {

    case ESP_A2D_CONNECTION_STATE_EVT: {
        bool is_connected =
            (param->conn_stat.state == ESP_A2D_CONNECTION_STATE_CONNECTED);
        ESP_LOGI(TAG, "%s (remote: %02x:%02x:%02x:%02x:%02x:%02x)",
                 is_connected ? "Connected" : "Disconnected",
                 param->conn_stat.remote_bda[0],
                 param->conn_stat.remote_bda[1],
                 param->conn_stat.remote_bda[2],
                 param->conn_stat.remote_bda[3],
                 param->conn_stat.remote_bda[4],
                 param->conn_stat.remote_bda[5]);

        if (is_connected && user_event_cb) {
            user_event_cb(BT_EVENT_CONNECTED);
        }

        if (!is_connected) {
            audio_running = false;
            if (user_event_cb) {
                user_event_cb(BT_EVENT_DISCONNECTED);
            }
            /* Re-enter discoverable mode so the next phone can connect. */
            esp_bt_gap_set_scan_mode(ESP_BT_CONNECTABLE,
                                     ESP_BT_GENERAL_DISCOVERABLE);
        }
        break;
    }

    case ESP_A2D_AUDIO_STATE_EVT:
        audio_running = (param->audio_stat.state == ESP_A2D_AUDIO_STATE_STARTED);
        ESP_LOGI(TAG, "Audio %s", audio_running ? "started" : "stopped");
        if (!audio_running) {
            xStreamBufferReset(pcm_stream);
        }
        break;

    case ESP_A2D_AUDIO_CFG_EVT: {
        esp_a2d_mcc_t *mcc = &param->audio_cfg.mcc;
        /* samp_freq is a bitmask, not a sequential enum.
         * The source picks a single rate during config,
         * so exactly one bit is set.
         * Priority: prefer higher rates. */
        uint8_t sf = mcc->cie.sbc_info.samp_freq;
        if      (sf & ESP_A2D_SBC_CIE_SF_48K) current_sample_rate = 48000;
        else if (sf & ESP_A2D_SBC_CIE_SF_44K) current_sample_rate = 44100;
        else if (sf & ESP_A2D_SBC_CIE_SF_32K) current_sample_rate = 32000;
        else if (sf & ESP_A2D_SBC_CIE_SF_16K) current_sample_rate = 16000;
        else                                  current_sample_rate = 44100;
        ESP_LOGI(TAG, "Audio cfg: SBC rate=%d Hz (samp_freq=0x%02x)",
                 current_sample_rate, sf);
        break;
    }

    default:
        break;
    }
}

/* ───────────────────────────────────────────────────────────────────
 *  A2DP data callback — runs in Bluedroid task context
 * ─────────────────────────────────────────────────────────────────── */

static void bt_a2dp_data_cb(const uint8_t *data, uint32_t len)
{
    if (!audio_running || len == 0) return;

    /* Push raw PCM bytes into the stream buffer. */
    size_t sent = xStreamBufferSend(pcm_stream, data, (size_t)len, 0);
    if (sent < len) {
        static uint32_t drop_count;
        drop_count++;
        if (drop_count < 10 || (drop_count % 100) == 0) {
            size_t avail = xStreamBufferBytesAvailable(pcm_stream);
            ESP_LOGW(TAG, "Drop #%lu: %u/%u bytes sent (buf avail=%u/%u)",
                     (unsigned long)drop_count,
                     (unsigned)len - (unsigned)sent, (unsigned)len,
                     (unsigned)avail, (unsigned)PCM_STREAM_BUF_SIZE);
        }

    }
}

/* ───────────────────────────────────────────────────────────────────
 *  DSP + I2S task — runs on the application core (core 1)
 * ─────────────────────────────────────────────────────────────────── */

static void dsp_i2s_task(void *arg)
{
    (void)arg;
    uint8_t buf[PCM_CHUNK];
    int rate = 44100;

    while (1) {
        /* Wait for audio data with a short timeout so we can
         * periodically check for BT event-driven SFX. */
        size_t got = xStreamBufferReceive(pcm_stream, buf, PCM_CHUNK,
                                          pdMS_TO_TICKS(50));

        /* Track sample rate changes negotiated by the BT source.
         * This must stay in sync with current_sample_rate so the
         * handler receives the correct rate for I2S/DSP init. */
        if (current_sample_rate != rate) {
            rate = current_sample_rate;
        }

        /* Let the user callback check for pending events (SFX, etc.)
         * even when no audio data arrived. */
        if (user_data_cb) {
            user_data_cb(got > 0 ? buf : NULL, (uint32_t)got, rate);
        }
    }
}
