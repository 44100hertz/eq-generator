/**
 * i2s_out.c — I2S DAC output driver for PCM5102A
 *
 * Uses ESP-IDF v5.x new I2S driver API.
 * The PCM5102A needs only BCK, LRCK (WS), and DIN — no MCLK.
 */

#include "i2s_out.h"
#include "driver/i2s_std.h"
#include "esp_log.h"
#include "freertos/FreeRTOS.h"

static const char *TAG = "i2s_out";

/* ── Pin assignment (generic ESP32 38-pin dev board) ──────────── */

#define I2S_BCK_PIN   26
#define I2S_LRCK_PIN  25
#define I2S_DATA_PIN  22

/* ── DMA buffer: ~20ms of stereo 16-bit at 48kHz ─────────────── */

#define I2S_DMA_DESC_NUM   4
#define I2S_DMA_FRAME_NUM  512

static i2s_chan_handle_t tx_handle = NULL;

/* ────────────────────────────────────────────────────────────────
 *  Init
 * ──────────────────────────────────────────────────────────────── */

int i2s_out_init(int sample_rate)
{
    esp_err_t err;

    /* Tear down previous instance if re-initializing at a new rate */
    if (tx_handle != NULL) {
        i2s_channel_disable(tx_handle);
        i2s_del_channel(tx_handle);
        tx_handle = NULL;
    }

    /* ── Allocate TX channel ─────────────────────────────────── */
    i2s_chan_config_t chan_cfg = {
        .id            = I2S_NUM_0,
        .role          = I2S_ROLE_MASTER,
        .dma_desc_num  = I2S_DMA_DESC_NUM,
        .dma_frame_num = I2S_DMA_FRAME_NUM,
        .auto_clear    = true,
    };

    err = i2s_new_channel(&chan_cfg, &tx_handle, NULL);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "i2s_new_channel failed: %d", err);
        return err;
    }

    /* ── Configure standard (Philips) mode ───────────────────── */
    i2s_std_config_t std_cfg = {
        .clk_cfg  = I2S_STD_CLK_DEFAULT_CONFIG(sample_rate),
        .slot_cfg = I2S_STD_PHILIPS_SLOT_DEFAULT_CONFIG(
                        I2S_DATA_BIT_WIDTH_16BIT,
                        I2S_SLOT_MODE_STEREO),
        .gpio_cfg = {
            .mclk = I2S_GPIO_UNUSED,
            .bclk = I2S_BCK_PIN,
            .ws   = I2S_LRCK_PIN,
            .dout = I2S_DATA_PIN,
            .din  = I2S_GPIO_UNUSED,
        },
    };

    err = i2s_channel_init_std_mode(tx_handle, &std_cfg);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "i2s_channel_init_std_mode failed: %d", err);
        i2s_del_channel(tx_handle);
        tx_handle = NULL;
        return err;
    }

    /* ── Enable the channel ──────────────────────────────────── */
    err = i2s_channel_enable(tx_handle);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "i2s_channel_enable failed: %d", err);
        i2s_del_channel(tx_handle);
        tx_handle = NULL;
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
    int16_t buf[2] = { left, right };
    size_t bytes_written = 0;

    esp_err_t err = i2s_channel_write(tx_handle, buf, sizeof(buf),
                                      &bytes_written, portMAX_DELAY);
    if (err != ESP_OK) {
        return err;
    }
    return 0;
}

int i2s_out_write(const int16_t *samples, uint32_t frame_count,
                  size_t *bytes_written)
{
    size_t nbytes = frame_count * 2 * sizeof(int16_t);
    size_t written = 0;
    esp_err_t err = i2s_channel_write(tx_handle, samples, nbytes,
                                      &written, portMAX_DELAY);
    if (bytes_written) *bytes_written = written;
    return (err == ESP_OK) ? 0 : (int)err;
}
