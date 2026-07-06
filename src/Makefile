# Build the harmonic bass enhancer shared library, LADSPA plugin,
# and native PipeWire filter.
#
# Usage:
#   make                  build enhancer.so (Python FFI)
#   make ladspa           build eqgen_ladspa.so (LADSPA plugin)
#   make filter           build eqgen_filter (native PipeWire filter)
#   make all              build everything
#   make install-ladspa   install LADSPA plugin to ~/.ladspa/
#   make clean            remove build artifacts
#   make test             build and run the C smoke test

CC       ?= gcc
CFLAGS   ?= -O2 -Wall -Wextra
LDFLAGS_SO ?= -shared -fPIC -lm

# ── PipeWire paths (NixOS) ────────────────────────────────────────
PW_DEV   := $(shell ls -d /nix/store/*pipewire*dev*/ 2>/dev/null | tail -1)
PW_RUN   := $(shell dirname $$(ldd /run/current-system/sw/bin/pw-cli 2>/dev/null | grep libpipewire | awk '{print $$3}') 2>/dev/null)
PW_CFLAGS := -I$(PW_DEV)include/pipewire-0.3 -I$(PW_DEV)include/spa-0.2
PW_LIBS   := -L$(PW_RUN) -lpipewire-0.3 -lm

SOURCES        = enhancer.c enhancer_api.c
HEADERS        = enhancer.h enhancer_api.h biquad_q28.h envelope.h dc_blocker.h
TARGET         = enhancer.so

LADSPA_SOURCES = enhancer.c ladspa/ladspa_wrapper.c
LADSPA_HEADERS = enhancer.h biquad_q28.h envelope.h dc_blocker.h eq_coeffs.h
LADSPA_TARGET  = eqgen_ladspa.so
LADSPA_INSTALL_DIR = $(HOME)/.ladspa

FILTER_SOURCES = enhancer.c ladspa/eqgen_filter.c
FILTER_HEADERS = enhancer.h biquad_q28.h envelope.h dc_blocker.h eq_coeffs.h
FILTER_TARGET  = eqgen_filter

TEST_SOURCES   = test_enhancer.c
TEST_TARGET    = test_enhancer

.PHONY: all ladspa filter install-ladspa clean test

all: $(TARGET) $(LADSPA_TARGET) $(FILTER_TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS_SO) -o $@ $(SOURCES)

$(LADSPA_TARGET): $(LADSPA_SOURCES) $(LADSPA_HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS_SO) -o $@ $(LADSPA_SOURCES)

ladspa: $(LADSPA_TARGET)

$(FILTER_TARGET): $(FILTER_SOURCES) $(FILTER_HEADERS)
	$(CC) $(CFLAGS) $(PW_CFLAGS) -o $@ $(FILTER_SOURCES) $(PW_LIBS)

filter: $(FILTER_TARGET)

install-ladspa: $(LADSPA_TARGET)
	mkdir -p $(LADSPA_INSTALL_DIR)
	cp $(LADSPA_TARGET) $(LADSPA_INSTALL_DIR)/

clean:
	rm -f $(TARGET) $(LADSPA_TARGET) $(FILTER_TARGET) $(TEST_TARGET)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

$(TEST_TARGET): $(TEST_SOURCES) $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DTEST_ENHANCER -o $@ $(TEST_SOURCES) $(SOURCES) -lm
