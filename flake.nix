{
  description = "MethylBin - metagenomic binning with neural networks";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            python311
            python311.pkgs.pip
            python311.pkgs.setuptools
            python311.pkgs.wheel

            cargo
            rustc
            linuxKernel.packages.linux_zen.perf
            cargo-cross
            cargo-release

            # For building
            openssl
            clang
            pkg-config
            libGL
            stdenv.cc.cc.lib
            glib

            bedtools
            prodigal
            hmmer
            samtools
          ];

          shellHook = ''
            export PYTHONPATH="$PWD:$PYTHONPATH"

            export LD_LIBRARY_PATH="${pkgs.libGL}/lib/:${pkgs.stdenv.cc.cc.lib}/lib/:${pkgs.glib.out}/lib/:$LD_LIBRARY_PATH"

            export MAMBA_ROOT_PREFIX="$PWD/.micromamba";
            if [ "$SHELL" = "/usr/bin/fish" ] || "$(basename "$SHELL" = "fish")"; then
              eval "$(micromamba shell hook --shell fish)"
            else
              eval "$(micromamba shell hook --shell bash)"
            fi
          '';
        };
        packages.default = pkgs.buildRustPackage {
          src = ./.;
          cargoLock = {
            lockFile = ./Cargo.lock;
          };
          buildInputs = [
            pkgs.openssl
            pkgs.pkg-config
          ];
          nativeBuildInputs = [
            pkgs.openssl_3_3
            pkgs.pkg-config
          ];
        };
      }
    );
}
