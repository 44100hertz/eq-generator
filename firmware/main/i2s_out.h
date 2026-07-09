/**
 * i2s_out.h — I2S DAC output driver for PCM5102A
 *
 * Configures ESP32 I2S0 in master mode for a standard external DAC.
 * Default pins work with generic 38-pin ESP32-WROOM dev boards.
 */

#pragma once
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Initialize I2S output.
 *
 *  Pins default to GPIO 26 (BCK), 25 (LRCK/WS), 22 (DATA).
 *  These are free on generic ESP32 dev boards and conflict with nothing common.
 *
 *  @param sample_rate  Usually 44100 or 48000, negotiated by A2DP.
 *  @return 0 on success, non-zero ESP-IDF error code on failure.
 */
int i2s_out_init(int sample_rate);

/** Write one stereo frame (2 × int16) to the I2S DMA buffer.
 *
 *  Non-blocking — data goes into a DMA ring buffer.
 *  Call from the A2DP data callback at audio rate.
 *
 *  @param left   Left channel sample, signed 16-bit.
 *  @param right  Right channel sample, signed 16-bit.
 *  @return 0 on success.
 */
int i2s_out_write_stereo(int16_t left, int16_t right);

/** Bulk-write multiple stereo frames to the I2S DMA buffer.
 *
 *  @param samples       Interleaved stereo int16_t buffer.
 *  @param frame_count   Number of stereo frames.
 *  @param bytes_written [out] actual bytes written (optional, may be NULL).
 *  @return 0 on success.
 */
int i2s_out_write(const int16_t *samples, uint32_t frame_count,
                  size_t *bytes_written);

#ifdef __cplusplus
}
#endif
