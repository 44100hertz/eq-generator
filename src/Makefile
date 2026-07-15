# Build the harmonic bass enhancer shared library and LADSPA plugin.
#
# Usage:
#   make                  build enhancer.so (Python FFI)
#   make ladspa           build eqgen_ladspa.so (LADSPA plugin)
#   make all              build everything
#   make install-ladspa   install LADSPA plugin to ~/.ladspa/
#   make clean            remove build artifacts
#   make test             build and run the C smoke test

CC       ?= gcc
CFLAGS   ?= -O2 -Wall -Wextra
LDFLAGS_SO ?= -shared -fPIC -lm

SOURCES        = enhancer.c enhancer_api.c
HEADERS        = enhancer.h enhancer_api.h biquad.h lp.h envelope.h eq_coeffs.h
TARGET         = enhancer.so

LADSPA_SOURCES = enhancer.c ladspa/ladspa_wrapper.c
LADSPA_HEADERS = enhancer.h biquad.h lp.h envelope.h eq_coeffs.h
LADSPA_TARGET  = eqgen_ladspa.so
LADSPA_INSTALL_DIR = $(HOME)/.ladspa

TEST_SOURCES   = test_enhancer.c
TEST_TARGET    = test_enhancer

FILTER_SOURCES = filter.c enhancer.c enhancer_api.c
FILTER_TARGET  = eqgen_filter

.PHONY: all ladspa install-ladspa filter clean test

all: $(TARGET) $(LADSPA_TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS_SO) -o $@ $(SOURCES)

$(LADSPA_TARGET): $(LADSPA_SOURCES) $(LADSPA_HEADERS)
	$(CC) $(CFLAGS) -I. $(LDFLAGS_SO) -o $@ $(LADSPA_SOURCES)

ladspa: $(LADSPA_TARGET)

install-ladspa: $(LADSPA_TARGET)
	mkdir -p $(LADSPA_INSTALL_DIR)
	cp $(LADSPA_TARGET) $(LADSPA_INSTALL_DIR)/

$(FILTER_TARGET): $(FILTER_SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -I. -o $@ $(FILTER_SOURCES) -lm

filter: $(FILTER_TARGET)

clean:
	rm -f $(TARGET) $(LADSPA_TARGET) $(TEST_TARGET) $(FILTER_TARGET)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

$(TEST_TARGET): $(TEST_SOURCES) $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DTEST_ENHANCER -o $@ $(TEST_SOURCES) $(SOURCES) -lm
