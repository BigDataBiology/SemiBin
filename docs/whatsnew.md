# What's New

## Unreleased (development version)
### User-visible improvements
- Users can now pass in the output of running mmseqs2 directly and SemiBin will
  use that instead of calling mmseqs itself (use option
  `--taxonomy-annotation-table`).
- The subcommand to generate cannot links is now called
  `generate_cannot_links`. The old name (`predict_taxonomy`) is kept as a
  deprecated alias.
- Provide pretrained models from soil, cat gut, human oral,pig gut, mouse gut,
  built environment, wastewater and global (training from all samples).
- Add `check_install` command and run `check_install` before easy command

### Bugfixes
- Fix bug with non-standard characters in sample names (#68).

## Version 0.5

*Released January 7 2022*

### User-visible improvements
- Reclustering is now the default (use `--no-recluster` to disable it; the
  option `--recluster` is deprecated and ignored) as the computational costs
  are much lower
- GTDB lazy downloading is now performed even if a non-standard directory is
  used
- The [CACHEDIR.TAG](https://bford.info/cachedir/) protocol was implemented
  (this is supported by several tools that perform tasks such as backups).

### Bugfixes
- Fix bug with `--min-len` (minimal length). Previously, only contigs greater
  than the given minimal length were used (instead of greater-equal to the
  minimal length).
- GTDB downloading was inconsistent in a few instances which have been fixed

### Internal improvements
- Much more efficient code (including lower memory usage) for binning,
  especially if a pretrained model is used. As an example, using a
  deeply-sequenced ocean sample, generating the data (`generate_data_single`
  step) goes down from 14 to 9 minutes; while binning (`bin` step, using
  `--recluster`) goes down from 10m17s (using 20GB of RAM, at peak) to 4m33
  (using 4.5 GB, at peak). Thus total time from BAM file to bins went down from
  25 to 14 minutes (using 4 threads) and peak RAM is now 4.5GB, making it
  usable on a typical laptop.

## Version 0.4.0

*Released 27 October 2021*

### User-visible improvements
- Add support for `.xz` FASTA files as input

### Internal improvements
- Removed BioPython dependency

### Bug fixes
- Fix bug when uncompressing FASTA files ([#42](https://github.com/BigDataBiology/SemiBin/issues/42))
- Fix bug when splitting data

## Version 0.3

*Released 10 August 2021*

### User-visible improvements
- Support training from several samples
- Remove `output_bin_path` if `output_bin_path` exists
- Make several internal parameters configuable: (1) minimum length of contigs to bin (`--min-len` parameter); (2) minimum length of contigs to break up in order to generate _must-link_ constraints (`--ml-threshold` parameter); (3) the ratio of the number of base pairs of contigs between 1000-2500 bp smaller than this value, the minimal length will be set as 1000bp, otherwise 2500bp (`--ratio` parameter).
- Add `-p` argument for `predict_taxonomy` mode

### Internal improvements
- Better code overall
- Fix `np.concatenate` warning
- Remove redundant matrix when clustering
- Better pretrained models
- Faster calculating dapth using Numpy
- Use correct number of threads in `kneighbors_graph()`

### Bugfixes

- Respect number of threads (`-p` argument) when training [(issue 34)](https://github.com/BigDataBiology/SemiBin/issues/34)

## Version 0.2

*Release 27 May 2021*

### User-visible improvements
- Change name to `SemiBin`
- Add support for training with several samples
- Test with Python 3.9
- Download mmseqs database with `--remove-tmp-file 1`
- Better output names
- Fix bugs when paths have spaces
- Fix installation issues by listing all the dependencies
- Add `download_GTDB` command
- Add `--recluster` option
- Add `--environment` option
- Add `--mode` option

### Internal improvements
- All around more robust code by including more error checking &amp; testing
- Better built-in models

## Version 0.1.1

*Released 21 March 2021*

**Bugfix release** fixing an issue with `minfasta-kbs`

## Version 0.1

*Released 21 March 2021*

- First release: testing version

