/**
 * sfx_player.c — Play embedded WAV sound effects through I2S.
 */

#include "sfx_player.h"
#include "rom/ets_sys.h"

void sfx_play(const int16_t *samples, uint32_t n,
              uint32_t rate_hz, uint32_t i2s_rate,
              int (*write_fn)(int16_t left, int16_t right))
{
    if (!samples || !n || !write_fn || i2s_rate == 0) return;

    /* Compute delay between samples at the I2S output rate.
     * Each iteration pushes one sample pair (stereo).  At e.g. 48000 Hz,
     * that's 48000 iterations per second → ~21 µs per iteration.
     *
     * We use a busy loop because the SFX are short (< 1 s) and
     * FreeRTOS task switching jitter would cause audible glitches. */
    uint32_t delay_us = (1000000UL + i2s_rate / 2) / i2s_rate;
    (void)rate_hz; /* documented only — playback uses i2s_rate */

    for (uint32_t i = 0; i < n; i++) {
        int16_t s = samples[i];
        write_fn(s, s);  /* mono → stereo */

        /* Spin-wait for one sample period.  ets_delay_us() is
         * an ESP-IDF ROM function available without includes. */
        if (delay_us > 1) {
            ets_delay_us(delay_us);
        }
    }
}
