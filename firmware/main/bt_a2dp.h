/**
 * bt_a2dp.h — Bluetooth A2DP sink
 *
 * Initializes Bluetooth Classic with Bluedroid, configures the ESP32
 * as an A2DP sink, and routes decoded PCM to a user callback.
 *
 * Usage:
 *   bt_a2dp_sink_init("eqgen", data_handler);
 *   // data_handler receives stereo int16 PCM at the negotiated rate.
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** PCM data callback type.
 *
 *  Called from a dedicated processing task at audio rate.
 *  Data is stereo interleaved int16, already decoded from SBC.
 *
 *  @param data   Pointer to raw PCM bytes (int16 interleaved).
 *  @param len    Byte count (multiple of 4 — two int16 per frame).
 *  @param rate   Negotiated sample rate (e.g. 44100 or 48000).
 */
typedef void (*bt_a2dp_data_cb_t)(const uint8_t *data, uint32_t len,
                                  int rate);

/** Initialize Bluetooth and start A2DP sink.
 *
 *  Device name appears during pairing.  Call once at boot.
 *  After this, the ESP32 is discoverable and will auto-accept
 *  incoming A2DP connections.
 *
 *  @param device_name  Name shown during Bluetooth pairing.
 *  @param data_cb      Called with decoded PCM frames.
 */
void bt_a2dp_sink_init(const char *device_name, bt_a2dp_data_cb_t data_cb);

#ifdef __cplusplus
}
#endif
