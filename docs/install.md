# Install

SemiBin can run on Python 3.6-3.9.

### Install from bioconda ###

```bash
conda create conda create -n SemiBin python==3.7
conda activate SemiBin
conda install -c bioconda semibin=0.2=pyh5e36f6f_1
```

### Install from source ###

Install dependence packages using conda: [MMseqs2](https://github.com/soedinglab/MMseqs2),[Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/),  [Fraggenescan](https://sourceforge.net/projects/fraggenescan/).

```bash
conda install -c conda-forge -c bioconda mmseqs2=13.45111
```
```bash
conda install -c bioconda bedtools hmmer fraggenescan
```

```bash
python setup.py install
```

