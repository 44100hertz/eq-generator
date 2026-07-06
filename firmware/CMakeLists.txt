# eqgen firmware — ESP32 Bluetooth A2DP sink with DSP
#
# Requires ESP-IDF v4.4+.  Build:
#   idf.py set-target esp32
#   idf.py build
#   idf.py flash
#
# Before building, generate coefficients:
#   python -m eqgen.cli.export -o firmware/main/eq_coeffs.h

cmake_minimum_required(VERSION 3.16)
include($ENV{IDF_PATH}/tools/cmake/project.cmake)
project(eqgen)
