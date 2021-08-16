# Install

SemiBin can run on Python 3.6-3.9.

### Install from bioconda ###

_Note_ : If you want to use SemiBin with GPU, you need to install Pytorch with GPU support. Or `conda install -c bioconda semibin` just install Pytorch with CPU support.

```bash
conda create conda create -n SemiBin python==3.7
conda activate SemiBin
conda install -c bioconda semibin
conda install pytorch torchvision torchaudio cudatoolkit=10.2 -c pytorch-lts
```

### Install from source ###

Install dependence packages using conda: [MMseqs2](https://github.com/soedinglab/MMseqs2),[Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/),  [Fraggenescan](https://sourceforge.net/projects/fraggenescan/).

```bash
conda install -c conda-forge -c bioconda mmseqs2=13.45111
```
```bash
conda install -c bioconda bedtools hmmer fraggenescan==1.30
```

```bash
python setup.py install
```

