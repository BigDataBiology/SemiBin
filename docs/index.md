# S<sup>3</sup>N<sup>2</sup>Bin

S<sup>3</sup>N<sup>2</sup>Bin is a command line tool for metagenomic binning with semi-supervised siamese neural network using additional information from reference genomes and contigs themselves. It will output the reconstructed bins in single sample/co-assembly/multi-samples binning mode.

### Single sample binning

Single sample binning means binning each sample into genomes after independently assembling. This mode allows for parallel binning of samples, but it do not use co-abundance information across samples.

### Co-assembly binning

Co-assembly binning means samples are co-assembled first and binning later.  This mode can generate better contigs and use co-abundance information, but co-assembly will lead to intersample chimeric contigs and binning based on co-assembly can not retain sample specific variation.

### Multi-samples binning

Multi-samples binning means multiple samples are assembled individually. Then these contigs are concatenated together and mapping reads to the concatenated contig from each sample. Then binning for every sample with abundance information from all samples. This mode can use  co-abundance information and retain the sample-specific variation at the same time. But this mode requires more computation resource. It takes more time when mapping reads back to the concatenated  fasta file.



## Commands

#### `single_easy_bin`

Reconstruct bins with single or co-assembly binning using one line command.

* `-i/--input-fasta` : Path to the input contig fasta file (gzip and bzip2 compression are accepted).
* `-b/--input-bam`: Path to the input BAM files.
* `-o/--output`: Output directory (will be created if non-existent).
* `--cannot-name:` Name for the cannot-link file(Default: cannot).
* `-r/--reference-db`: GTDB reference file.(Default: $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB) If not set `--reference-db` and can not find GTDB in the default path, we will download GTDB to the default path.
* `-p/--processes/-t/--threads`: Number of CPUs used(0: use whole).
* `--minfasta-kbs`: minimum bin size in kilo-basepairs (Default: 200).
* `--epoches`: Number of epoches used in the training process(Default: 20).
* `--batch-size`: Batch size used in the training process(Default: 2048).
* `--max-node`: Percentage of contigs that considered to be binned(Default: 1).
* `--max-edges`: The maximum number of edges that can be connected to one contig(Default: 200).

#### multi_easy_bin

Reconstruct bins with multi-samples binning using one line command.


* `-b/--input-bam`: Path to the input BAM files. Unlike in `single_easy_bin`, you can pass multiple BAM files, one per sample
* `-s/--separator`: Used when multiple samples binning to separate sample name and contig name(Default is `:`).

The following options (including synonyms) are the same as for
`single_easy_bin`: `--input-fasta`, `--output`, `--reference-db`,
`--processes`, `--minfasta-kbs`, `--epoches`, `--batch-size`, `--max-node`, and
`--max-edges`.



#### predict_taxonomy

Run the contig annotations using mmseqs with GTDB and generate cannot-link file used in the semi-supervised deep learning model training.

* `-i/--input-fasta` : Path to the input contig fasta file (gzip and bzip2 compression are accepted).
* `-o/--output`: Output directory (will be created if non-existent).
* `--cannot-name:` Name for the cannot-link file(Default: cannot).
* `-r/--reference-db`: GTDB reference file.(Default: $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB) If not set `--reference-db` and can not find GTDB in the default path, we will download GTDB to the default path.

#### generate_data_single

Generate training data(data.csv;data_split.csv) for single and co-assembly binning.

* `-i/--input-fasta` : Path to the input contig fasta file (gzip and bzip2 compression are accepted).
* `-b/--input-bam`: Path to the input bam files. If there are several samples, you can input several bam files.
* `-o/--output`: Output directory (will be created if non-existent).
* `-p/--processes/-t/--threads`: Number of CPUs used(0: use whole).

#### generate_data_multi

Generate training data(data.csv;data_split.csv) for multi-samples binning.

* `-i/--input-fasta` : Path to the input contig fasta file (gzip and bzip2 compression are accepted).
* `-b/--input-bam`: Path to the input bam files. If there are several samples, you can input several bam files.
* `-o/--output`: Output directory (will be created if non-existent).
* `-s/--separator`: Used when multiple samples binning to separate sample name and contig name(Default is `:`).
* `-p/--processes/-t/--threads`: Number of CPUs used(0: use whole).

#### bin

Training the model and clustering contigs into bins.

* `-i/--input-fasta` : Path to the input contig fasta file (gzip and bzip2 compression are accepted).
* `-b/--input-bam`: Path to the input bam files. If there are several samples, you can input several bam files.
* `-o/--output`: Output directory (will be created if non-existent).
* `--data`: Path to the input data.csv file.
* `--data_split`: Path to the input data_split.csv file.
* `-c/--cannot-link` : Path to the input cannot link file generated from other additional biological information, one row for each cannot link constraint. The file format: contig_1,contig_2.
* `--minfasta-kbs`: minimum bin size in kilo-basepairs (Default: 200).
* `--epoches`: Number of epoches used in the training process(Default: 20).
* `--batch-size`: Batch size used in the training process(Default: 2048).
* `--max-node`: Percentage of contigs that considered to be binned(Default: 1).
* `--max-edges`: The maximum number of edges that can be connected to one contig(Default: 200).
* `-p/--processes/-t/--threads`: Number of CPUs used(0: use whole).



