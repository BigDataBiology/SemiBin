# S<sup>3</sup>N<sup>2</sup>Bin (Semi-supervised Siamese Neural Network for metagenomic binning)


_NOTE_: This tool is still in development. You are welcome to try it out and
feedback is appreciated, but expect some bugs/rapid changes until it
stabilizes.

Command tool for metagenomic binning with semi-supervised deep learning using
information from reference genomes.

## Install

S<sup>3</sup>N<sup>2</sup>Bin runs on Python 3.6-3.8 and requires [Bedtools](https://github.com/arq5x/bedtools2) to be available for calculating depth from bam files.

### Install from source

You can download the source code from github and install with the standard

```bash
python setup.py install
```

## Generate Cannot-link contrains

You can use [mmseqs](https://github.com/soedinglab/MMseqs2) or [CAT](https://github.com/dutilh/CAT) (or other contig annotation tools) to get taxonomic classifications of contigs. Then you can use the script `script/concatenate.py` to generate the cannot-link file(contig1, contig2)
that can be used in S<sup>3</sup>N<sup>2</sup>Bin.

#### mmseqs

```bash
python script/concatenate.py -i taxonomy.tsv -c contig.fna -s sample-name -o output --mmseqs
```

#### CAT

```bash
python script/concatenate.py -i CAT.out -c contig.fna -s sample-name -o output --CAT
```

## Examples

### Easy-bin mode 

You can get the results with one line code. Easy-bin command can be used in single-sample and co-assembly binning modes(contig annotations using mmseqs with GTDB reference genome).

```bash
S3N2Bin easy-bin -i contig.fna -b *.bam --GTDB-path /mmseqs_data/GTDB -o output
```

If you do not set the path of GTDB, S<sup>3</sup>N<sup>2</sup>Bin will download GTDB  to your output folder.

### Advance-bin mode(Need to generate cannot-link file before running)

#### Single sample/co-assembly binning

```bash
S3N2Bin advance-bin -i contig.fna -b *.bam -c cannot-link.txt -o output 
```

#### Multiple samples binning (Must set `-s` parameter)

```bash
S3N2Bin advance-bin -i whole_contig.fna -b *.bam -c *.txt -s C -o output
```

#### Multi-samples binning pipeline

(1) Concatenate all contigs from all samples together. Make sure that names of samples are unique and id for every contig is <sample_name><\separator><contig_id>. Concatenated contig format is:

```bash
<S1>O<Contig_1>
ATGCAAAA
<S1>O<Contig_2>
ATGCAAAA
<S1>O<Contig_3>
ATGCAAAA
<S2>O<Contig_1>
ATGCAAAA
<S2>O<Contig_2>
ATGCAAAA
<S3>O<Contig_1>
ATGCAAAA
<S3>O<Contig_2>
ATGCAAAA
<S3>O<Contig_3>
ATGCAAAA
```

(2) Map reads to the concatenated contig file to get the bam files.

(3) Generate cannot-link files for every sample. The name of the cannot-link file should be <sample_name>.txt. Make sure the sample name here is same to that in step 1.

(4) Run S<sup>3</sup>N<sup>2</sup>Bin.

For more details(i.e. --split-run, which is reconmended for projects with large samples), [read the docs](https://s3n2bin.readthedocs.io/en/latest/usage/). 

## Output

The output folder will contain

1. Datasets used for training and clustering.

2. Saved semi-supervised deep learning model.

3. Output bins.

4. Some intermediate files.

When single sample/co-assembly binning, reconstructed bins are in `output_recluster_bins` directory. When multi-samples binning, reconstructed bins from all samples are in `bins` directory. 

For more details, [read the docs](https://s3n2bin.readthedocs.io/en/latest/output/). 

