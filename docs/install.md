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


```bash
conda install -c conda-forge -c bioconda mmseqs2=13.45111
```
```bash
conda install -c bioconda bedtools hmmer fraggenescan
```

```bash
python setup.py install
```

