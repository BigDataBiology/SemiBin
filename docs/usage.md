# Usage

## Easy single/co-assembly binning mode

You can get the results with one line of code. The `single_easy_bin` command can be used in
single-sample and co-assembly binning modes (contig annotations using mmseqs
with GTDB reference genome). `single_easy_bin` includes the following parts:
`predict_taxonomy`,`generate_data_single` and `bin`.

Your inputs consist of a (1) a contig FASTA file and (2) a BAM file of mapped
reads to the contigs (see [Generating inputs to S³N²Bin](generate.html)).

Run S<sup>3</sup>N<sup>2</sup>Bin in `single_easy_bin` mode.

```bash
S3N2Bin single_easy_bin -i contig.fna -b *.bam -r /mmseqs_data/GTDB -o output
```

If you do not set the path of GTDB, S<sup>3</sup>N<sup>2</sup>Bin will download GTDB  to $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB.

## Easy multi-samples binning mode

The `multi_easy_bin` command can be used in
multi-samples binning modes (contig annotations using mmseqs
with GTDB reference genome). `multi_easy_bin` includes following parts: `predict_taxonomy`, `generate_data_multi` and `bin`.

(1) Concatenate all contigs from all samples together. Make sure that names of
samples are unique and rename the contigs to `<sample_name>:<contig_id>`
(_i.e._, the sample name, `:`, and then the within-sample contig name; if your
sample names contain a colon, you can use the `--separator` argument to use a
different separator). For example, your concatenated FASTA file
(`contig_whole.fna`) could look like the following

```bash
>S1:Contig_1
AGATAATAAAGATAATAATA
>S1:Contig_2
CGAATTTATCTCAAGAACAAGAAAA
>S1:Contig_3
AAAAAGAGAAAATTCAGAATTAGCCAATAAAATA
>S2:Contig_1
AATGATATAATACTTAATA
>S2:Contig_2
AAAATATTAAAGAAATAATGAAAGAAA
>S3:Contig_1
ATAAAGACGATAAAATAATAAAAGCCAAATCCGACAAAGAAAGAACGG
>S3:Contig_2
AATATTTTAGAGAAAGACATAAACAATAAGAAAAGTATT
>S3:Contig_3
CAAATACGAATGATTCTTTATTAGATTATCTTAATAAGAATATC
```

(2) Map reads to the `contig_whole.fna`.

(3) Run S<sup>3</sup>N<sup>2</sup>Bin with `multi_easy_bin` mode.

```bash
S3N2Bin multi_easy_bin -i contig_whole.fna -b *.bam -o output
```

See above for the comment about the location where the GTDB database is stored.

## Advanced-bin mode

Especially for multi-samples binning, running with `multi_easy_bin` takes much time when processing samples serially. You can split the step and manually run S³N²Bin parallel.

### Generate cannot-link constrains

S³N²Bin has builtin support for
[mmseqs2](https://github.com/soedinglab/MMseqs2) (the default) and
[CAT](https://github.com/dutilh/CAT) for generating contig taxonomic
classifications and generating cannot-link file. See below for format
specifications if you want to use another tool.

#### Contig annotation with mmseqs (GTDB reference):

To use mmseqs (default in S<sup>3</sup>N<sup>2</sup>Bin and producing the best
results in benchmarks), you can use subcommand `predict_taxonomy` to generate
the cannot-link file. If you want to use CAT, you have run CAT first anduse
the script `script/concatenate.py` to generate the cannot-link file(contig1,
contig2) that can be used in S<sup>3</sup>N<sup>2</sup>Bin.

```bash
S3N2Bin predict_taxonomy -i contig.fna -o output
```

Default GTDB path is $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB. If you do not set
`--reference-db` and do not find GTDB in the default path. S³N²Bin will
download GTDB to the default path.

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

### Generating data for training and clustering (`data.csv;data_split.csv`)

#### Single/co-assembly binning

```bash
S3N2Bin generate_data_single -i contig.fna -b *.bam -o output
```

#### Multi-samples binning

```bash
S3N2Bin generate_data_multi -i contig_whole.fna -b *.bam -o output
```

### Binning (training and clustering)

```bash
S3N2Bin bin -i contig.fna --data data.csv --data-split data_split.csv -c cannot.txt -o output
```

If a GPU is available, S<sup>3</sup>N<sup>2</sup>Bin will automatically take
advantage of it.

