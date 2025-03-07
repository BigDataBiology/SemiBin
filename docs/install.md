# Install

SemiBin can run on Python 3.7-3.13.

## Install with pixi

The current recommended way to install SemiBin is to use [pixi](https://pixi.sh/). Pixi will use the packages from conda-forge and bioconda to install SemiBin and its dependencies.

### CPU-only mode

Create a `pixi.toml` file with the following content.

```toml
[project]
authors = ["Luis Pedro Coelho <luis@luispedro.org>"]
channels = ["conda-forge", "bioconda"]
name = "semibin_install"
platforms = ["linux-64", "osx-64"]
version = "0.1.0"

[tasks]

[dependencies]
semibin = ">=2.1.0,<3"
```

Then, run `pixi install` in the same directory as the `pixi.toml` file to download and install SemiBin2.


### With GPU support

If you want to use SemiBin with GPU, you need to install Pytorch with GPU support as well. Starting with the example above, you need to add `pytorch-gpu` to the `dependencies` section and `cuda` to the `system-requirements` section.

```toml
[project]
authors = ["Luis Pedro Coelho <luis@luispedro.org>"]
channels = ["conda-forge", "bioconda"]
name = "semibin_install"
platforms = ["linux-64"]
version = "0.1.0"

[tasks]

[dependencies]
semibin = ">=2.1.0,<3"
pytorch-gpu = "*"

[system-requirements]
cuda = "12.0"
```

## Install from bioconda


Pixi is now the recommended way to install SemiBin. However, if you prefer to use conda, you can install SemiBin with it

### Simple Mode

```bash
conda create -n SemiBin
conda activate SemiBin
conda install -c conda-forge -c bioconda semibin
```

### GPU mode from conda

If you want to use SemiBin with GPU, you need to install Pytorch with GPU support.

```bash
conda create -n SemiBin
conda activate SemiBin
conda install -c conda-forge -c bioconda semibin
conda install pytorch torchvision torchaudio cudatoolkit=10.2 -c pytorch-lts
```

## Install from source

You will need the following dependencies:
- [Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/)
- [Prodigal](https://github.com/hyattpd/Prodigal)


You can obtain them from conda with the following commands

```bash
conda install -c bioconda bedtools hmmer samtools
```

Then, installing should be a simple matter of running:

```bash
pip install .
```

## Alternative ways of running SemiBin

If you use one of these pipelines, we ask that you cite both the pipeline author and the [SemiBin manuscript](https://www.nature.com/articles/s41467-022-29843-y) (as well as any other pipeline-embedded tools which contribute to your results).

- [ATLAS](https://metagenome-atlas.github.io/) is a Snakemake-based pipeline for metagenomics, which includes SemiBin (as well as other binners and tools).
- [Galaxy toolshed](https://toolshed.g2.bx.psu.edu/view/iuc/suite_semibin/) also includes a [SemiBin wrapper](https://toolshed.g2.bx.psu.edu/view/iuc/suite_semibin/)


