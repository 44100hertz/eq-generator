/**
 * sfx_player.h — Play embedded WAV sound effects through I2S.
 *
 * Generated sfx_data.h contains int16 PCM arrays (e.g. sfx_connected[]).
 * Call sfx_play() from the dsp_i2s task context — it blocks until done.
 */

#pragma once
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Play a PCM array through I2S.
 *
 *  @param samples   Pointer to int16 mono sample array.
 *  @param n         Number of samples.
 *  @param rate_hz   Sample rate the array was recorded at (for timing).
 *  @param i2s_rate  Current I2S output rate (samples are played as-is at this rate).
 *  @param write_fn  Callback: void write_fn(int16_t left, int16_t right)
 *                   Called once per sample pair (mono → stereo duplicate).
 */
void sfx_play(const int16_t *samples, uint32_t n,
              uint32_t rate_hz, uint32_t i2s_rate,
              int (*write_fn)(int16_t left, int16_t right));

#ifdef __cplusplus
}
#endif
