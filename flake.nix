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

            micromamba

            bedtools
            hmmer
            samtools
          ];

          shellHook = ''
            export PYTHONPATH="$PWD:$PYTHONPATH"

            export MAMBA_ROOT_PREFIX="$PWD/.micromamba";
            if [ "$SHELL" = "/usr/bin/fish" ] || "$(basename "$SHELL" = "fish")"; then
              eval "$(micromamba shell hook --shell fish)"
            else
              eval "$(micromamba shell hook --shell bash)"
            fi
          '';
        };
      }
    );
}
