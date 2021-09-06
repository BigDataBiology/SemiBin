# SemiBin

SemiBin is a command line tool for metagenomic binning with semi-supervised siamese neural network using additional information from reference genomes and contigs themselves. It will output the reconstructed bins in single sample/co-assembly/multi-samples binning mode.

### Single sample binning

Single sample binning means that each sample is assembled and binned
independently.

This mode allows for parallel binning of samples and avoid cross-sample
chimeras, but it does not use co-abundance information across samples.

### Co-assembly binning

Co-assembly binning means samples are co-assembled first (as if the pool of
samples were a single sample) and binned later.

This mode can generate better contigs (especially from species that are at a
low abundance in any individual sample) and use co-abundance information, but
co-assembly can lead to intersample chimeric contigs and binning based on
co-assembly dows not retain sample specific variation. It is appropriate when
the samples are very similar.

### Multi-sample binning

With multi-sample binning, multiple samples are assembled and binned
individually, but _information from multiple samples is used together_.
This mode can use co-abundance information and retain sample-specific
variation at the same time. However, it has increased computational costs.

This mode is implemented by concatenating the contigs assembled from the
individual samples together and then mapping reads from each sample to this
concatenated database.

## Commands

#### `single_easy_bin`

Reconstruct bins with single or co-assembly binning using one line command.

* `-i/--input-fasta` : Path to the input contig fasta file (gzip and bzip2 compression are accepted).
* `-b/--input-bam`: Path to the input BAM files. You can pass multiple BAM files, one per sample.
* `-o/--output`: Output directory (will be created if non-existent).
* `--cannot-name:` Name for the cannot-link file (Default: cannot).
* `-r/--reference-db`: GTDB reference file.(Default: $HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB) If not set `--reference-db` and can not find GTDB in the default path, we will download GTDB to the default path.
* `-p/--processes/-t/--threads`: Number of CPUs used(0: use whole).
* `--minfasta-kbs`: minimum bin size in kilo-basepairs (Default: 200).
* `--recluster` : Recluster bins(which will take more time and better results).
* `--epoches`: Number of epoches used in the training process(Default: 20).
* `--batch-size`: Batch size used in the training process(Default: 2048).
* `--max-node`: Percentage of contigs that considered to be binned(Default: 1).
* `--max-edges`: The maximum number of edges that can be connected to one contig(Default: 200).
* `--random-seed`: Random seed to reproduce results.
* `--environment`: Environment for the built-in model(human_gut/dog_gut/ocean).
* `--ratio` : If the ratio of the number of base pairs of contigs between 1000-2500 bp smaller than this value, the minimal length will be set as 1000bp, otherwise2500bp. If you set -m parameter, you do not need to use this parameter. If you use SemiBin with multi steps and you use this parameter, please use this parameter consistently with all subcommands(Default: 0.05). 
* `-m/--min-len` : Minimal length for contigs in binning. If you use SemiBin with multi steps and you use this parameter, please use this parameter consistently with all subcommands.(Default: SemiBin chooses 1000bp or 2500bp according the ratio of the number of base pairs of contigs between 1000-2500 bp).
* `--ml-threshold` : Length threshold for generating must-link constraints.(By default, the threshold is calculated from the contig, and the default minimum value is 4,000 bp)

#### multi_easy_bin

Reconstruct bins with multi-samples binning using one line command.

* `-b/--input-bam`: Path to the input BAM files. You can pass multiple BAM files, one per sample.
* `-s/--separator`: Used when multiple samples binning to separate sample name and contig name(Default is `:`).

The following options (including synonyms) are the same as for
`single_easy_bin`: `--input-fasta`, `--output`, `--reference-db`,
`--processes`, `--minfasta-kbs`, `--recluster`,`--epoches`, `--batch-size`, `--max-node`, and
`--max-edges`, `--random-seed`, `--ratio`, `--min-len`, `--ml-threshold`.

#### predict_taxonomy

Run the contig annotations using mmseqs with GTDB and generate cannot-link file used in the semi-supervised deep learning model training.

The following options are the same as for `single_easy_bin`: `-i/--input-fasta`, `-o/--output`, `--cannot-name`, `-r/--reference-db`, `--ratio`, `--min-len` and `--ml-threshold`.

#### generate_data_single

Generate training data(data.csv;data_split.csv) for single and co-assembly binning.

The following options are the same as for `single_easy_bin`: `-i/--input-fasta`,  `-b/--input-bam`, `-o/--output`, `-p/--processes/-t/--threads`, `--ratio`, `--min-len`, `--ml-threshold`.

#### generate_data_multi

Generate training data(data.csv;data_split.csv) for multi-samples binning.

The following options are the same as for `single_easy_bin`: `-i/--input-fasta`, `-o/--output`, `-p/--processes/-t/--threads`, `--ratio`, `--min-len`, `--ml-threshold`.

The following options are the same as for `multi_easy_bin`: `-b/--input-bam`, `-s/--separator`.

#### train ####

Training the model.

* `--data`: Path to the input data.csv file.
* `--data_split`: Path to the input data_split.csv file.
* `-c/--cannot-link` : Path to the input cannot link file generated from other additional biological information, one row for each cannot link constraint. The file format: contig_1,contig_2.
* `--mode`:  [single/several] Train models from one sample or several samples(train model across several samples can get better pre-trained model for single-sample binning.) In several mode, must input data, data_split, cannot, fasta files for corresponding sample with same order. *Note:* You can just set `several` with this option when single-sample binning. Training from several samples with multi-sample binning is not support.

The following options are the same as for `single_easy_bin`: `-i/--input-fasta`,  `-o/--output`, `--epoches`, `--batch-size`, `-p/--processes/-t/--threads`, `--random-seed`, `--ratio`, `--min-len`.

#### bin

Clustering contigs into bins.

* `--model`: Path to the trained model.

The following options are the same as for `single_easy_bin`: `-i/--input-fasta`, `-o/--output`, `--minfasta-kbs`, `--recluster`, `--max-node`, `--max-edges`, `-p/--processes/-t/--threads`, `--random-seed`, `--environment`, `--ratio`, `--min-len`.

#### download_GTDB

Download reference genomes(GTDB).

* `-r/--reference-db`: GTDB reference file path to download(~/path).(Default: $HOME/.cache/SemiBin/mmseqs2-GTDB) If not set `--reference-db` , we will download GTDB to the default path.

