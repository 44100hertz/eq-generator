{
  description = "Speaker EQ correction suite for JamesDSP (Python + EEL)";

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
            echo "EQ correction suite ready."
            echo "  Generate EQ:   python eqgen.py -m <meas.wav> -t <target.wav> [-n <noise.wav>] -o <output>"
            echo "  Run tests:     ./run_analysis.sh"
            echo "  List tests:    ./run_analysis.sh --list"
            echo "  Sanity check:  ./sanitycheck.sh"
          '';
        };
      });
}
