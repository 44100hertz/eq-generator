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

#include <string.h>
#include "freertos/FreeRTOS.h"
#include "freertos/stream_buffer.h"
#include "freertos/task.h"

#include "esp_log.h"
#include "nvs_flash.h"
#include "esp_bt.h"
#include "esp_bt_main.h"
#include "esp_bt_device.h"
#include "esp_a2dp_api.h"
#include "esp_avrc_api.h"

static const char *TAG = "bt_a2dp";

/* ── StreamBuffer between BT task and DSP task ────────────────────── */

/* Hold ~100ms of stereo 16-bit audio at 48 kHz (worst case).
 * 100ms × 48000 × 4 bytes/frame ≈ 19200 → round up to 20480.     */
#define PCM_STREAM_BUF_SIZE  20480

static StreamBufferHandle_t pcm_stream;
static bt_a2dp_data_cb_t    user_data_cb;
static int                  current_sample_rate = 44100;
static bool                 audio_running = false;

/* ── Forward declarations ──────────────────────────────────────────── */

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

    /* ── A2DP sink ────────────────────────────────────────────── */
    ESP_ERROR_CHECK(esp_a2d_register_callback(bt_a2dp_event_cb));
    ESP_ERROR_CHECK(esp_a2d_sink_register_data_callback(bt_a2dp_data_cb));
    ESP_ERROR_CHECK(esp_a2d_sink_init());

    /* ── StreamBuffer + DSP task (runs on app core) ───────────── */
    pcm_stream = xStreamBufferCreate(PCM_STREAM_BUF_SIZE, 4);
    xTaskCreatePinnedToCore(dsp_i2s_task, "dsp_i2s", 4096, NULL, 5,
                            NULL, 1);

    /* ── Set discoverable ─────────────────────────────────────── */
    esp_bt_gap_set_scan_mode(ESP_BT_CONNECTABLE,
                             ESP_BT_GENERAL_DISCOVERABLE);

    ESP_LOGI(TAG, "Ready as \"%s\". Waiting for connection…", device_name);
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
            (param->conn_state.state == ESP_A2D_CONNECTION_STATE_CONNECTED);
        ESP_LOGI(TAG, "%s (remote: %02x:%02x:%02x:%02x:%02x:%02x)",
                 is_connected ? "Connected" : "Disconnected",
                 param->conn_state.remote_bda[0],
                 param->conn_state.remote_bda[1],
                 param->conn_state.remote_bda[2],
                 param->conn_state.remote_bda[3],
                 param->conn_state.remote_bda[4],
                 param->conn_state.remote_bda[5]);

        if (!is_connected) {
            audio_running = false;
            /* Re-enter discoverable mode so the next phone can connect. */
            esp_bt_gap_set_scan_mode(ESP_BT_CONNECTABLE,
                                     ESP_BT_GENERAL_DISCOVERABLE);
        }
        break;
    }

    case ESP_A2D_AUDIO_STATE_EVT:
        audio_running = (param->audio_state.state == ESP_A2D_AUDIO_STATE_STARTED);
        ESP_LOGI(TAG, "Audio %s", audio_running ? "started" : "stopped");
        if (!audio_running) {
            xStreamBufferReset(pcm_stream);
        }
        break;

    case ESP_A2D_AUDIO_CFG_EVT: {
        esp_a2d_audio_cfg_t *a = &param->audio_cfg;
        ESP_LOGI(TAG, "Codec=%s rate=%u ch=%u bits=%u",
                 (a->mcc.type == ESP_A2D_MCT_SBC) ? "SBC" : "?",
                 (unsigned)a->mcc.cie.sbc.sample_freq,
                 (unsigned)a->mcc.cie.sbc.ch_mode,
                 (unsigned)a->mcc.cie.sbc.alloc);

        /* Map SBC sample_freq enum to Hz.
         * 0=16k, 1=32k, 2=44.1k, 3=48k */
        switch (a->mcc.cie.sbc.sample_freq) {
            case 0:  current_sample_rate = 16000; break;
            case 1:  current_sample_rate = 32000; break;
            case 2:  current_sample_rate = 44100; break;
            case 3:  current_sample_rate = 48000; break;
            default: current_sample_rate = 44100; break;
        }
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

    /* Push raw PCM bytes into the stream buffer (non-blocking).
     * If the DSP task is falling behind, we skip this frame — no
     * point buffering old audio. The stream buffer is sized to
     * handle normal scheduling jitter.                                  */
    size_t sent = xStreamBufferSend(pcm_stream, data, (size_t)len, 0);
    if (sent < len) {
        /* DSP task is overloaded — drop silently.  We'll catch up
         * once the buffer drains.                                       */
        static uint32_t drop_count;
        if (++drop_count < 10 || (drop_count % 100) == 0) {
            ESP_LOGW(TAG, "Drop: %u/%u bytes (total drops: %lu)",
                     (unsigned)len - (unsigned)sent, (unsigned)len,
                     drop_count);
        }
        /* Reset the buffer to recover. */
        xStreamBufferReset(pcm_stream);
    }
}

/* ───────────────────────────────────────────────────────────────────
 *  DSP + I2S task — runs on the application core (core 1)
 * ─────────────────────────────────────────────────────────────────── */

#define PCM_CHUNK 512   /* bytes to read per tick (128 stereo frames) */

static void dsp_i2s_task(void *arg)
{
    (void)arg;
    uint8_t buf[PCM_CHUNK];
    int rate = 44100;

    while (1) {
        /* Wait for at least one stereo frame (4 bytes).
         * On first connect, or after a disconnect, block until
         * audio arrives.  If audio is already streaming, this
         * returns immediately.                                            */
        size_t got = xStreamBufferReceive(pcm_stream, buf, PCM_CHUNK,
                                          portMAX_DELAY);

        /* Re-init I2S if the sample rate changed. */
        if (current_sample_rate != rate) {
            rate = current_sample_rate;
            ESP_LOGI(TAG, "Sample rate changed to %d — re-initializing I2S", rate);
            /* Re-init I2S handled by calling i2s_out_init again.
             * We call it via the user callback's first invocation
             * with the new rate.                                          */
        }

        if (user_data_cb) {
            user_data_cb(buf, (uint32_t)got, rate);
        }
    }
}
