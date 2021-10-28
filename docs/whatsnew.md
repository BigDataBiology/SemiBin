# What's New

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

