# Main component — links DSP sources from ../../src/ and our firmware code.
#
# Before building, generate coefficients from measurements:
#   python -m eqgen.cli.wire build <speaker>    ← uses measurements/<speaker>/
# See firmware/README.md for details.

set(COEFFS_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../../src/eq_coeffs.h")
if(NOT EXISTS "${COEFFS_HEADER}")
    message(FATAL_ERROR
        "src/eq_coeffs.h not found.\n"
        "Run: python -m eqgen.cli.wire build <speaker>\n"
        "See firmware/README.md for details.")
endif()

idf_component_register(
    SRCS
        "main.c"
        "i2s_out.c"
        "bt_a2dp.c"
        "sfx_player.c"
        "../../src/enhancer.c"
    INCLUDE_DIRS
        "."
        "../../src"
    REQUIRES
        nvs_flash
        driver
        esp_adc
        esp_driver_i2s
        bt
)
