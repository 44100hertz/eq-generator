# Main component — links DSP sources from ../../src/ and our firmware code.

idf_component_register(
    SRCS
        "main.c"
        "i2s_out.c"
        "bt_a2dp.c"
        "../../src/enhancer.c"
    INCLUDE_DIRS
        "."
        "../../src"
    REQUIRES
        nvs_flash
        driver
        esp_adc
)
