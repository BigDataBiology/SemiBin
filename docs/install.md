# Install

SemiBin can run on Python 3.7-3.10.

## Install from bioconda

### Simple Mode

The simplest way to install is to use [conda](https://conda.io/).

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
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- (optionally) [Fraggenescan](https://sourceforge.net/projects/fraggenescan/)


You can obtain them from conda with the following commands

```bash
conda install -c conda-forge -c bioconda mmseqs2=13.45111
conda install -c bioconda bedtools hmmer fraggenescan samtools
```

Then, installing should be a simple matter of running:

```bash
python setup.py install
```

## Alternative ways of running SemiBin

If you use one of these pipelines, we ask that you cite both the pipeline author and the [SemiBin manuscript](https://www.nature.com/articles/s41467-022-29843-y) (as well as any other pipeline-embedded tools which contribute to your results).

- [ATLAS](https://metagenome-atlas.github.io/) is a Snakemake-based pipeline for metagenomics, which includes SemiBin (as well as other binners and tools).
- [Galaxy toolshed](https://toolshed.g2.bx.psu.edu/view/iuc/suite_semibin/) also includes a [SemiBin wrapper](https://toolshed.g2.bx.psu.edu/view/iuc/suite_semibin/)


