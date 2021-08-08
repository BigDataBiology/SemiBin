# What's New

## Version 0.3

*Release 10 August 2021*

### User-visible improvements
- Remove output_bin_path if output_bin_path exists
- Support train from several samples
- Add --min-len
- Add --ratio
- Add --ml-threshold
- Add -p for SemiBin predict_taxonomy

### Internal improvements
- Better code
- Fix np.concatenate warning
- Remove redundant matrix when clustering
- Better pretrained models
- Faster calculating dapth using Numpy
- Add  -p in kneighbors_graph()

### Bugfixes

- Fix bug -p does not work when training [(Issue 34)](https://github.com/BigDataBiology/SemiBin/issues/34)

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

