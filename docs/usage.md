# Usage

## Easy single/co-assembly binning mode

You can get the results with one line of code. The `single_easy_bin` command can be used in
single-sample and co-assembly binning modes (contig annotations using mmseqs
with GTDB reference genome). `single_easy_bin` includes the following parts:
`predict_taxonomy`,`generate_data_single` and `bin`.

Your inputs consist of a (1) a contig FASTA file and (2) a BAM file of mapped
reads to the contigs (see [Generating inputs to SemiBin](generate.html)).

Run SemiBin in `single_easy_bin` mode.

```bash
SemiBin single_easy_bin -i contig.fna -b *.bam -r /mmseqs_data/GTDB -o output 
```

If you do not set the path of GTDB, SemiBin will download GTDB  to $HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB. You can set `--recluster` to use the reclustering part with single-copy genes described in the paper.

You can use `--environment` with(human_gut, dog_gut and ocean) to use our built-in model. (**Note:** Recommended way, which will save much time for contig annotations and model training, and also get very good results) 

```bash
SemiBin single_easy_bin -i contig.fna -b *.bam -o output --environment human_gut
```

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

(3) Run SemiBin with `multi_easy_bin` mode.

```bash
SemiBin multi_easy_bin -i contig_whole.fna -b *.bam -o output
```

See above for the comment about the location where the GTDB database is stored. You can set `--recluster` to use the reclustering part with single-copy genes described in the paper.

## Advanced-bin mode

Especially for multi-samples binning, running with `multi_easy_bin` takes much time when processing samples serially. You can split the step and manually run SemiBin parallelly. 

Another advantage is that you do not annotate contigs and train models for every sample. You can use our built-in trained models and get a good model from your own datasets for single-sample binning. Then you can transfer this to your datasets to get binning results and it can save much time.

### Generate cannot-link constrains

SemiBin has built-in support for
[MMseqs2](https://github.com/soedinglab/MMseqs2) (the default) and
[CAT](https://github.com/dutilh/CAT) for generating contig taxonomic
classifications and generating cannot-link file. See below for format
specifications if you want to use another tool.

#### Contig annotation with MMseqs (GTDB reference):

To use MMseqs (default in SemiBin and producing the best
results in benchmarks), you can use subcommand `predict_taxonomy` to generate
the cannot-link file. If you want to use CAT, you have run CAT first anduse
the script `script/concatenate.py` to generate the cannot-link file(contig1,
contig2) that can be used in SemiBin.

```bash
SemiBin predict_taxonomy -i contig.fna -o output
```

Default GTDB path is $HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB. If you do not set
`--reference-db` and do not find GTDB in the default path. SemiBin will
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
SemiBin generate_data_single -i contig.fna -b *.bam -o output
```

#### Multi-samples binning

```bash
SemiBin generate_data_multi -i contig_whole.fna -b *.bam -o output
```

### Training

```bash
SemiBin train -i contig.fna --data data.csv --data-split data_split.csv -c cannot.txt -o output --mode single
```

### Binning

If you want to use our provided models(human_gut/dog_gut/ocean) and you have a model from your samples, you do not need to generate cannot-link constrains and train model, you can just generate data.csv and get the binning results. 

```bash
SemiBin bin -i contig.fna --model model.h5 --data data.csv -o output 
```

or

```bash
SemiBin bin -i contig.fna --data data.csv -o output --environment human_gut 
```

#### Getting a model from your project with single-sample binning

For example, you can subsample several samples(i.e. 5) as training samples and several samples as testing samples. Then you can train models from every training sample and test the models in the testing samples. Finally you can use the best model in other samples and get the binning results.

(1) Generate data.csv/data_split.csv for every sample

```bash
SemiBin generate_data_single -i contig.fna -b *.bam -o output
```

(2) Generate cannot-link for every sample

```bash
SemiBin predict_taxonomy -i contig.fna -o output
```

(3) Train a pre-trained model across several samples

```bash
SemiBin train -i *.fna --data *.csv --data-split *.csv -c cannot*.txt -o output --mode several
```

(4) Bin with the trained model

```bash
SemiBin bin -i contig.fna --model model.h5 --data data.csv -o output 
```

If a GPU is available, SemiBin will automatically take
advantage of it.

### Command(whole)

```bash
SemiBin predict_taxonomy -i contig.fna -o output
```

```bash
SemiBin generate_data_single -i contig.fna -b *.bam -o output
SemiBin generate_data_multi -i contig_whole.fna -b *.bam -o output
```

```bash
SemiBin train -i contig.fna --data data.csv --data-split data_split.csv -c cannot.txt -o output --mode single
```

```bash
SemiBin bin -i contig.fna --model model.h5 --data data.csv -o output 
```

### Command(with pre-trained model) ###

```bash
SemiBin generate_data_single -i contig.fna -b *.bam -o output
SemiBin generate_data_multi -i contig_whole.fna -b *.bam -o output
```

```bash
SemiBin bin -i contig.fna --model model.h5 --data data.csv -o output 
SemiBin bin -i contig.fna --data data.csv -o output 
```



