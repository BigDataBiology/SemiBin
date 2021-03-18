# Usage

## Easy single/co-assembly binning mode

You can get the results with one line of code. The `single_easy_bin` command can be used in
single-sample and co-assembly binning modes (contig annotations using mmseqs
with GTDB reference genome). `single_easy_bin` includes the following parts: `predict_taxonomy`,`generate_data_single` and `bin`.

(1) Mapping reads to the contig fasta file. 

Example using Bowtie2 (but another read mapping tool would also work):

```bash
bowtie2-build -f contig.fna contig.fna -p 16

bowtie2 -q --fr contig.fna -1 reads_1.fq.gz -2 reads_2.fq.gz -S contig.sam -p 64

samtools view -h -b -S contig.sam -o contig.bam -@ 64

samtools view -b -F 4 contig.bam -o contig.mapped.bam -@ 64

samtools sort -m 1000000000 contig.mapped.bam -o contig.mapped.sorted.bam -@ 64

samtools index contig.mapped.sorted.bam
```

(2) Run S<sup>3</sup>N<sup>2</sup>Bin with single_easy_bin mode.

```bash
S3N2Bin single_easy_bin -i contig.fna -b *.bam -r /mmseqs_data/GTDB -o output
```

If you do not set the path of GTDB, S<sup>3</sup>N<sup>2</sup>Bin will download GTDB  to $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB.

## Easy multi-samples binning mode

The `multi_easy_bin` command can be used in
multi-samples binning modes (contig annotations using mmseqs
with GTDB reference genome). `multi_easy_bin` includes following parts: `predict_taxonomy`, `generate_data_multi` and `bin`.

(1) Concatenate all contigs from all samples together. Make sure that names of samples are unique and id for every contig is <sample_name>:<contig_id>. ':' is the separator that used to split the contig name to sample_name and contig_id. You can use any separator you want by set `--separator`(Default is  `:`). Just make sure that the separator will not introduce confusion when splitting. Concatenated contig format is:

```bash
<S1>:<Contig_1>
ATGCAAAA
<S1>:<Contig_2>
ATGCAAAA
<S1>:<Contig_3>
ATGCAAAA
<S2>:<Contig_1>
ATGCAAAA
<S2>:<Contig_2>
ATGCAAAA
<S3>:<Contig_1>
ATGCAAAA
<S3>:<Contig_2>
ATGCAAAA
<S3>:<Contig_3>
ATGCAAAA
```

(2) Mapping reads to the contig_whole.fna . 

(3) Run S<sup>3</sup>N<sup>2</sup>Bin with multi_easy_bin mode.

```bash
S3N2Bin multi_easy_bin -i contig_whole.fna -b *.bam -r /mmseqs_data/GTDB -o output -s :
```

## Advanced-bin mode

Especially for multi-samples binning, running with `multi_easy_bin` takes much time when processing samples serially. You can split the step and manually run S³N²Bin parallel.

### Generate cannot-link constrains

S³N²Bin has builtin support for
[mmseqs2](https://github.com/soedinglab/MMseqs2) (the default) and
[CAT](https://github.com/dutilh/CAT) for generating contig taxonomic classifications and generating cannot-link file. See below for format specifications if you want to use another tool.

#### Contig annotation with mmseqs (GTDB reference):

Default GTDB path is $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB. If you do not set `--reference-db` and do not find GTDB in the default path. S³N²Bin will download GTDB to the default path.

```bash
S3N2Bin predict_taxonomy -i contig.fna -r /mmseqs_data/GTDB -o output
```

#### Contig annotation with CAT

```bash
CAT contigs \
        -c contig.fasta \
        -d CAT_prepare_20200304/2020-03-04_CAT_database \
        --path_to_prodigal path_to_prodigal \
        --path_to_diamond path_to_diamond \
        -t CAT_prepare_20200304/2020-03-04_taxonomy \
        -o CAT_output/CAT \
        --force \
        -f 0.5 \
        --top 11 \
        --I_know_what_Im_doing \
        --index_chunks 1

CAT add_names \
    CAT_output/CAT.contig2classification.txt \
    -o CAT_output/CAT.out \
    -t CAT_prepare_20200304/2020-03-04_taxonomy \
    --force \
    --only_official
```

Generate cannot-link constrains

```bash
python script/concatenate.py -i CAT.out -c contig.fasta -s sample-name -o output --CAT
```

### Generating data for training and clustering(data.csv;data_split.csv)

#### Single/co-assembly binning

```bash
S3N2Bin generate_data_single -i contig.fna -b *.bam -o output
```

#### Multi-samples binning

```bash
S3N2Bin generate_data_multi -i contig_whole.fna -b *.bam -o output -s :
```

### Binning(training and clustering)

If you run S<sup>3</sup>N<sup>2</sup>Bin on a GPU server, S<sup>3</sup>N<sup>2</sup>Bin will run on GPU automatically.

```bash
S3N2Bin bin -i contig.fna --data data.csv --data-split data_split.csv -c cannot.txt -o output
```