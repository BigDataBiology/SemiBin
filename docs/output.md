# Outputs of SemiBin

## Single sample/co-assembly binning

* `output_bins`: directory of the final reconstructed bins (after reclustering by default).
* `model.h5`: saved deep learning model (self-supervised by default).
* `data.csv/data_split.csv`: data used in the training of deep learning model.
* `*_data_cov.csv/*_data_split_cov.csv`: coverage data generated from depth file.
* `cannot/cannot.txt`: cannot-link file used in training (only present when using semi-supervised mode).
* `bins_info.tsv`: table with basic information on each bin (name, total number of basepairs, number of contigs, N50, and L50; more columns may be added in the future).
* `contig_bins.tsv`: table mapping each contig to its assigned bin.

If `--write-pre-reclustering-bins` is passed, pre-reclustering bins are also written to a separate `output_prerecluster_bins` directory.

## Multi-samples binning

* `bins`: Reconstructed bins from all samples.
* `samples/*.fasta`: Contig fasta file for every sample from the input whole_contig.fna.
* `samples/*_data_cov.csv`: same as in single sample/coassembly binning.
* `samples/{sample-name}/`: directory of the output of SemiBin for every sample (same as that in single sample/coassembly binning).

