{
  description = "autoeq + EEL_VM CLI toolchain";

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
        packages = {
          eel-cli = pkgs.stdenv.mkDerivation {
            pname = "eel-cli";
            version = "0.1.0";

            src = ./EEL_VM;

            nativeBuildInputs = with pkgs; [ gcc ];

            buildPhase = ''
              cd ..

              # All .c source files from the VS project
              SRCS="
                cpthread.c
                nseel-compiler.c
                nseel-ram.c
                loose_eel.c
                fft.c
                s_str.c
                y.tab.c
                numericSys/codelet.c
                numericSys/cpoly.c
                numericSys/FFTConvolver.c
                numericSys/FilterDesign/cos_fib_paraunitary.c
                numericSys/FilterDesign/eqnerror.c
                numericSys/FilterDesign/firls.c
                numericSys/FilterDesign/generalFdesign.c
                numericSys/FilterDesign/polyphaseASRC.c
                numericSys/FilterDesign/polyphaseFilterbank.c
                numericSys/HPFloat/atox.c
                numericSys/HPFloat/constant.c
                numericSys/HPFloat/cxaop.c
                numericSys/HPFloat/cxbasic.c
                numericSys/HPFloat/cxconstant.c
                numericSys/HPFloat/cxconvf.c
                numericSys/HPFloat/cxexp.c
                numericSys/HPFloat/cxhypb.c
                numericSys/HPFloat/cxidiv.c
                numericSys/HPFloat/cxpow.c
                numericSys/HPFloat/cxprcmp.c
                numericSys/HPFloat/cxtrig.c
                numericSys/HPFloat/hpaconf.c
                numericSys/HPFloat/prcxpr.c
                numericSys/HPFloat/print.c
                numericSys/HPFloat/prxpr.c
                numericSys/HPFloat/sfmod.c
                numericSys/HPFloat/shift.c
                numericSys/HPFloat/xadd.c
                numericSys/HPFloat/xchcof.c
                numericSys/HPFloat/xdiv.c
                numericSys/HPFloat/xevtch.c
                numericSys/HPFloat/xexp.c
                numericSys/HPFloat/xfmod.c
                numericSys/HPFloat/xfrac.c
                numericSys/HPFloat/xhypb.c
                numericSys/HPFloat/xivhypb.c
                numericSys/HPFloat/xivtrg.c
                numericSys/HPFloat/xlog.c
                numericSys/HPFloat/xmul.c
                numericSys/HPFloat/xneg.c
                numericSys/HPFloat/xprcmp.c
                numericSys/HPFloat/xpwr.c
                numericSys/HPFloat/xsigerr.c
                numericSys/HPFloat/xsqrt.c
                numericSys/HPFloat/xtodbl.c
                numericSys/HPFloat/xtoflt.c
                numericSys/HPFloat/xtrig.c
                numericSys/libsamplerate/samplerate.c
                numericSys/libsamplerate/src_linear.c
                numericSys/libsamplerate/src_sinc.c
                numericSys/MersenneTwister.c
                numericSys/quadprog.c
                numericSys/SolveLinearSystem/inv.c
                numericSys/SolveLinearSystem/mldivide.c
                numericSys/SolveLinearSystem/mrdivide.c
                numericSys/SolveLinearSystem/pinv.c
                numericSys/SolveLinearSystem/qr_fact.c
                numericSys/solvopt.c
              "

              gcc -O2 -I. -lm $SRCS -o eel_CLI
            '';

            installPhase = ''
              mkdir -p $out/bin
              cp eel_CLI $out/bin/
            '';
          };
        };

        defaultPackage = self.packages.${system}.eel-cli;

        devShells.default = pkgs.mkShell {
          buildInputs = [ self.packages.${system}.eel-cli ];
          shellHook = ''
            echo "EEL_VM CLI ready — run 'eel_CLI <script.eel>' to test scripts"
          '';
        };
      });
}
