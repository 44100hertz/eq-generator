# eqgen -- Speaker EQ Correction Suite
#
# Usage:
#   make                  build enhancer.so (Python FFI)
#   make all              build enhancer.so + LADSPA
#   make test             run all Python tests
#   make graph-check      full pipeline report → HTML graphs (via ARGS)
#   make clean            remove build artifacts
#   make eqgen            design EQ curve (JSON)
#   make audition         audition through C DSP -> WAV
#   make wire-setup       live PipeWire system-wide EQ
#   make wire-teardown    remove PipeWire wiring
#   make export           export coeffs for ESP32
#   make server           start web UI server (preset management)

.PHONY: all clean test graph-check eqgen audition wire-setup wire-teardown export server flash

# -- C DSP build -------------------------------------------------------

all:
	$(MAKE) -C src all

clean:
	$(MAKE) -C src clean
	rm -rf output/

# -- Python tests ------------------------------------------------------

test:
	python -m eqgen.tests.run_all

graph-check:
	python -m eqgen.cli.graph_check $(ARGS)

# -- Pipeline entry points ---------------------------------------------
# Usage: make eqgen ARGS="-m meas.wav -t target.wav -o eq.json"
eqgen:
	python -m eqgen.cli.eqgen $(ARGS)

audition:
	python -m eqgen.cli.audition $(ARGS)

wire-setup:
	python -m eqgen.cli.wire setup $(ARGS)

wire-teardown:
	python -m eqgen.cli.wire teardown

export:
	python -m eqgen.cli.export $(ARGS)

server:
	@PID_FILE=output/server.pid; \
	if [ -f "$$PID_FILE" ]; then \
		kill $$(cat "$$PID_FILE") 2>/dev/null || true; \
		rm -f "$$PID_FILE"; \
		sleep 0.5; \
	fi
	python -m eqgen.server $(ARGS)

# -- ESP32 firmware ----------------------------------------------------
# Usage: make flash ARGS="technics/standing"
flash:
	python -m eqgen.cli.wire build $(ARGS)
	cd firmware && idf.py build flash
