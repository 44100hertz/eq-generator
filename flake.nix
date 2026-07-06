{
  description = "Speaker EQ correction suite for JamesDSP (Python + C DSP)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            (pkgs.python3.withPackages (ps: with ps; [
              numpy
              scipy
              matplotlib
            ]))
          ];
          shellHook = ''
            echo "eqgen — speaker EQ correction suite"
            echo "  Design EQ:   python -m eqgen.cli.eqgen -m <meas.wav> -t <target.wav> -o eq.json"
            echo "  Audition:    python -m eqgen.cli.audition <speaker> /tmp/out --tracks song.flac"
            echo "  Live wire:   python -m eqgen.cli.wire setup <speaker>"
            echo "  Export:      python -m eqgen.cli.export --speaker small -o src/eq_coeffs.h"
            echo "  Run tests:   make test"
            echo "  Graph check: make graph-check"
          '';
        };
      });
}
