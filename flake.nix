{
  description = "Speaker EQ correction suite (Python + C DSP + ESP32 firmware)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs-esp-dev.url = "github:mirrexagon/nixpkgs-esp-dev";
  };

  outputs = { self, nixpkgs, flake-utils, nixpkgs-esp-dev }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        esp-idf-esp32 = nixpkgs-esp-dev.packages.${system}.esp-idf-esp32;
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            (pkgs.python3.withPackages (ps: with ps; [
              numpy
              scipy
              matplotlib
              fastapi
              uvicorn
              evdev
            ]))
            pkgs.gdb
            pkgs.libnotify
            esp-idf-esp32
          ];
          shellHook = ''
            echo "eqgen — speaker EQ correction suite"
            echo "  Design EQ:   python -m eqgen.cli.eqgen -m <meas.wav> -t <target.wav> -o eq.json"
            echo "  Audition:    python -m eqgen.cli.audition <speaker> /tmp/out --tracks song.flac"
            echo "  Live wire:   python -m eqgen.cli.wire setup <speaker>"
            echo "  Export:      python -m eqgen.cli.export --speaker small -o src/eq_coeffs.h"
            echo "  ESP32 build: cd firmware && idf.py build"
            echo "  ESP32 flash: cd firmware && idf.py flash"
            echo "  Run tests:   make test"
            echo "  Graph check: make graph-check"
          '';
        };
      });
}
