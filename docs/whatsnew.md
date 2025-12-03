# What's New

## Unreleased changes

### User-visible changes

- Add support for `SEMIBIN_DEBUG` environment variable to enable debug logging (overrides command line flags)

## Version 2.2.0

*Released Mar 20, 2025*

This is a maintenance release with many small improvement rather than a single big new feature. Upgrading is recommended, but not crucial.

### User-visible changes
- Remove `SemiBin` command. Only `SemiBin1` and `SemiBin2` are available (and `SemiBin1` is deprecated). The only reason to use `SemiBin1` is if you have old scripts that use it. It will be removed in the next release.
- Better logging: Always log to file in DEBUG level and log command-line arguments. Print version number in logs.
- Better error messages in several instances
- check_install: Prints out information on the GPU

### Deprecations
- SemiBin: Deprecate `--prodigal-output-faa` argument
- No longer check for `mmseqs` in `check_install` (it is not a hard requirement)

### Internal improvements and bugfixes
- Respect the number of threads requested better ([#140](https://github.com/BigDataBiology/SemiBin/issues/140))
- SemiBin: Better method to save the model which is more compatible with newer versions of PyTorch. Added a subcommand to update old models to the new format (`update_model`)
- SemiBin: Switch to pixi for testing (and recommend it in the README/[installation](install) instructions)
- Convert to `pyproject.toml` instead of `setup.py`
- Do not fail if no bins are produced ([#170](https://github.com/BigDataBiology/SemiBin/issues/170) &amp; [#173](https://github.com/BigDataBiology/SemiBin/issues/173))

## Version 2.1.0

*Released Mar 6, 2024*

Main new feature is adding support for using output of strobealign-aemb.

Use of the `SemiBin` command (instead of `SemiBin2`) will continue to work, but
print a warning and set a delay to ask users to upgrade.

### User-visible changes

- Support running SemiBin with [strobealign-aemb](https://github.com/ksahlin/strobealign/releases/tag/v0.13.0) (`--abundance`/`-a`)
- Add `citation` subcommand
- Introduce separate `SemiBin1` command as use of `SemiBin` is now deprecated and will trigger a warning

### Internal improvements
- Code simplification and refactor
- deprecation: Deprecate --orf-finder=fraggenescan option
- Update abundance normalization

### Bugfixes
- SemiBin: do not use more processes than can be taken advantage of [#155](https://github.com/BigDataBiology/SemiBin/issues/155)

## Version 2.0.2

*Released Oct 31, 2023*

### Bugfix release

Fixes issue with `multi_easy_bin --write-pre-reclustering-bins` [#128 on GH](https://github.com/BigDataBiology/SemiBin/issues/128#issuecomment-1784742272)

## Version 2.0.1

*Released Oct 21, 2023*

This is a bugfix release for _version 2.0.0_.

## Version 2.0.0

*Released Oct 20, 2023*

### User-visible changes

- Running SemiBin now writes a log file in the output directory
- The `concatenate_fasta` subcommand now supports compression
- Adds `bin_short` subcommand as alias for `bin` (by analogy with `bin_long`)


## Version 1.5.1 (SemiBin2 beta)

*Released Mar 7, 2023*

### Bugfixes

- Fix use of `--no-recluster` with multi_easy_bin ([#128](https://github.com/BigDataBiology/SemiBin/issues/128)).

## Version 1.5.0 (SemiBin2 beta)

*Released Jan 17, 2023*

Big change is the addition of a `SemiBin2` script, which is still experimental, but should be a slightly nicer interface.
See [[upgrading to SemiBin2](semibin2)]

### User-visible improvements

- Added a new option for ORF finding, called `fast-naive` which is an internal very fast implementation.
- Added the possibility of bypassing ORF finding altogether by providing prodigal outputs directly (or any other gene prediction in the right format)
- Command line argument checking is more exhaustive instead of exiting at first error
- Added `--quiet` flag to reduce the amount of output printed
- Better `--help` (group required arguments separately)
- Add `--output-compression` option to compress outputs
- Add `--tag-output` option which allows for control of the output filenames (and also makes the anvi'o compatible — see discussion at [#123](https://github.com/BigDataBiology/SemiBin/issues/123).
- Add contig->bin mapping table ([#123](https://github.com/BigDataBiology/SemiBin/issues/123))
- `SemiBin.main.main1` and `SemiBin.main.main2` can now be called as a function with command line arguments (`main1` corresponds to _SemiBin1_ and `main2` corresponds to _SemiBin2_)

```python
import SemiBin.main

...

SemiBin.main.main2(['single_easy_bin', '--input-fasta', ...])
```

## Version 1.4.0: long reads binning!

*Released December 15, 2022*

Big change is the added binning algorithm for assemblies from long-read datasets.

The overall structure of the pipeline is still similar to what was [manuscript](https://www.nature.com/articles/s41467-022-29843-y), but when clustering, it does not use infomap, but another procedure (an iterative version of DBSCAN).

Use the flag `--sequencing-type=long_read` to enable an alternative clustering that works better with long reads.

### Other user-visible improvements

- Better error checking at multiple steps in the pipeline so that processes that will crash are caught as early as possible
- Add `--allow-missing-mmseqs2` flag to `check_install` subcommand (eventually, self-supervision will be the default and mmseqs2 will be an optional dependency)

### Command line parameter deprecations

The previous arguments should continue to work, but going forward, the newer arguments are probably a better API.

- Selecting self-supervised learning is now done with the `--self-supervised` flag (instead of `--training-type=self`)
- Training from multiple samples is now enabled with the `--train-from-many` flag (instead of `--mode=several`)

### Bugfixes

- The output table sometimes had the wrong path in `v1.3`. This has been fixed
- Prodigal is now run in a more robust manner when using multiple threads ([#106](https://github.com/BigDataBiology/SemiBin/issues/106))

## Version 1.3.1

*Release December 9, 2022*

### Bugfixes

- Made `--training-type` argument optional (defaults to `semi` to keep backwards compatibility)


## Version 1.3.0

*Released November 4 2022*

### User visible improvements

- Added _self-supervised learning mode_ (see [[Training SemiBin models](training)] for more details)

### Bugfixes

- Fix output table to contain correct paths
- Fix mispelling in argument name `--epochs` (the old variation, `--epoches` is still accepted for backwards compatibility, but should be considered deprecated)

## Version 1.2.0

*Released October 19 2022*

### User visible improvements

- Pretrained model from chicken caecum (contributed by [Florian Plaza Oñate](https://scholar.google.com/citations?hl=zh-CN&user=-gE5y_4AAAAJ&view_op=list_works&sortby=pubdate))
- Output table with basic information on bins (including N50 & L50)
- When reclustering is used (default), output the unreclusted bins into a directory called `output_prerecluster_bins`
- Added `--verbose` flag and silenced some of the output when it is not used
- Use coloredlogs (if package is available)

## Version 1.1.1

*Released September 27 2022*

### Bugfixes

- Completely remove use of `atomicwrites` package ([#97](https://github.com/BigDataBiology/SemiBin/issues/97))

## Version 1.1.0

*Released September 21 2022*

### User-visible improvements

- Support .cram format input ([#104](https://github.com/BigDataBiology/SemiBin/issues/104))
- Support using depth file from Metabat2 ([#103](https://github.com/BigDataBiology/SemiBin/issues/103))
- More flexible specification of prebuilt models (case insensitive, normalize `-` and `_`)
- Better output message when no bins are produced

### Bugfixes

- Fix bug using `atomicwrite` on certain network filesystems ([#97](https://github.com/BigDataBiology/SemiBin/issues/97))

### Internal improvements

- Remove torch version restriction (and test on Python 3.10)


## Version 1.0.3

*Released August 3 2022*

### Bugfixes

- Fix coverage parsing when value is not an integer ([#103](https://github.com/BigDataBiology/SemiBin/issues/103))
- Fix multi_easy_bin with taxonomy file given on the command line (see discussion at [#102](https://github.com/BigDataBiology/SemiBin/issues/102))


## Version 1.0.2

*Released July 8 2022*

### Bugfixes

- Fix ([#93](https://github.com/BigDataBiology/SemiBin/issues/93)) more
  thoroughly ([#101](https://github.com/BigDataBiology/SemiBin/issues/101))

## Version 1.0.1

*Released May 9 2022*

### Bugfixes

- Fix edge case when calling prodigal with more threads than contigs
  ([#93](https://github.com/BigDataBiology/SemiBin/issues/93))

## Version 1.0.0

*Released April 29 2022*

This coincides with the publication of the
[manuscript](https://www.nature.com/articles/s41467-022-29843-y).

### User-visible improvements

- More balanced file split when calling prodigal in parallel should take better advantage of multiple threads
- Fix bug when long stretches of Ns are present ([#87](https://github.com/BigDataBiology/SemiBin/issues/87)]
- Better error messages
  ([#90](https://github.com/BigDataBiology/SemiBin/issues/90) &amp;
   [#91](https://github.com/BigDataBiology/SemiBin/issues/91)])

### Bugfixes

- Fix bugs in training from multiple samples
- Fix bug in incorporating CAT results

## Version 0.7

*Released March 2 2022*

This release solves [issues running on Mac OS X](https://github.com/BigDataBiology/SemiBin/issues/77).

### User-visible improvements

- Improved `check_install` command: it now prints out paths and correctly handles optionality of FragGeneScan/prodigal
- Add `concatenate_fasta` command to combine fasta files for multi-sample binning
- Add option `--tmpdir` to set temporary directory
- Substitute FragGeneScan with Prodigal (FragGeneScan can still be used with `--orf-finder` parameter). FragGeneScan caused issues, especially on Mac OSX

### Internal improvements
- Reuse `markers.hmmout` file to make the training from several samples faster

## Version 0.6

*Released February 7 2022*

### User-visible improvements
- Provide pretrained models from soil, cat gut, human oral,pig gut, mouse gut,
  built environment, wastewater and global (training from all samples).
- Users can now pass in the output of running mmseqs2 directly and SemiBin will
  use that instead of calling mmseqs itself (use option
  `--taxonomy-annotation-table`).
- The subcommand to generate cannot links is now called
  `generate_cannot_links`. The old name (`predict_taxonomy`) is kept as a
  deprecated alias.
- Similarly, sequence features (_k_-mer and abundance) are generated using the
  commands `generate_sequence_features_single` and
  `generate_sequence_features_multi` (for single- and multi-sample modes,
  respectively). The old names (`generate_data_single`/`generate_data_multi`)
  are kept as deprecated aliases.
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

