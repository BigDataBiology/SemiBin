# Install

S<sup>3</sup>N<sup>2</sup>Bin can run on Python 3.6-3.8.

Install from source

Install dependence packages using conda: [Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/),  [Fraggenescan](https://sourceforge.net/projects/fraggenescan/) and [cmake](https://cmake.org/).

```bash
conda install -c bioconda bedtools hmmer fraggenescan
```

```bash
conda install -c anaconda cmake=3.19.6
```

```bash
python setup.py install
```
