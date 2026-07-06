/**
 * i2s_out.c — I2S DAC output driver for PCM5102A
 *
 * Uses ESP-IDF's standard I2S driver (compatible with v4.4+).
 * The PCM5102A needs only BCK, LRCK, and DIN — no MCLK.
 */

#include "i2s_out.h"
#include "driver/i2s.h"
#include "esp_log.h"

static const char *TAG = "i2s_out";

/* ── Pin assignment (generic ESP32 38-pin dev board) ──────────── */

#define I2S_BCK_PIN   26
#define I2S_LRCK_PIN  25
#define I2S_DATA_PIN  22

/* ── DMA buffer: holds ~20ms of stereo 16-bit audio at 48kHz ────
 *  2 buffers × 1024 samples × 2 channels × 2 bytes = 4096 bytes each.
 *  Safe margin for Bluetooth task scheduling jitter.
 */

#define I2S_DMA_BUF_COUNT  4
#define I2S_DMA_BUF_LEN    1024

/* ────────────────────────────────────────────────────────────────
 *  Init
 * ──────────────────────────────────────────────────────────────── */

int i2s_out_init(int sample_rate)
{
    i2s_config_t cfg = {
        .mode                 = I2S_MODE_MASTER | I2S_MODE_TX,
        .sample_rate          = (uint32_t)sample_rate,
        .bits_per_sample      = I2S_BITS_PER_SAMPLE_16BIT,
        .channel_format       = I2S_CHANNEL_FMT_RIGHT_LEFT,
        .communication_format = I2S_COMM_FORMAT_STAND_I2S,
        .intr_alloc_flags     = ESP_INTR_FLAG_LEVEL1,
        .dma_buf_count        = I2S_DMA_BUF_COUNT,
        .dma_buf_len          = I2S_DMA_BUF_LEN,
        .use_apll             = false,
        .tx_desc_auto_clear   = true,
        .fixed_mclk           = 0,
    };

    i2s_pin_config_t pins = {
        .bck_io_num           = I2S_BCK_PIN,
        .ws_io_num            = I2S_LRCK_PIN,
        .data_out_num         = I2S_DATA_PIN,
        .data_in_num          = I2S_PIN_NO_CHANGE,
    };

    esp_err_t err;

    err = i2s_driver_install(I2S_NUM_0, &cfg, 0, NULL);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "i2s_driver_install failed: %d", err);
        return err;
    }

    err = i2s_set_pin(I2S_NUM_0, &pins);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "i2s_set_pin failed: %d", err);
        i2s_driver_uninstall(I2S_NUM_0);
        return err;
    }

    /* PCM5102A starts outputting on the first BCK edge — no mute/unmute needed. */
    ESP_LOGI(TAG, "I2S0 master TX @ %d Hz on BCK=%d LRCK=%d DATA=%d",
             sample_rate, I2S_BCK_PIN, I2S_LRCK_PIN, I2S_DATA_PIN);

    return 0;
}

/* ────────────────────────────────────────────────────────────────
 *  Write one stereo frame
 * ──────────────────────────────────────────────────────────────── */

int i2s_out_write_stereo(int16_t left, int16_t right)
{
    /* Pack two int16 samples into a single 32-bit word.
     * I2S transmits MSB first: right channel goes to the high half-word. */
    int16_t buf[2] = { left, right };
    size_t bytes_written = 0;

    esp_err_t err = i2s_write(I2S_NUM_0, buf, sizeof(buf),
                              &bytes_written, portMAX_DELAY);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "i2s_write failed: %d", err);
        return err;
    }
    return 0;
}
