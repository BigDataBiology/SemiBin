# SemiBin2

Starting with version 1.5 (officially _SemiBin2 beta_), installing the SemiBin package installs two scripts: `SemiBin` and `SemiBin2`.
They have the same functionality, but slightly different interfaces.
As of version 2.0.0, the older `SemiBin` command is _not recommended_ (except for backwards compability) and newer projects should use `SemiBin2`.

## Upgrading to SemiBin2

1. If you are using the `easy_*` workflows, then they will probably continue to
   work exactly the same (except that you will get better results faster).
2. Outputs are now **always** in a directory called `output_bins` (unless you explicitly ask for the pre-reclustered bins to be written out with the `--write-pre-reclustering-bins` option).
3. By default, bins are in file named as `SemiBin_{label}.fa.gz` (and
   compressed with _gzip_ as the name indicates; you can change the compression with the `--compression` flag).

Points `2` and `3` may require some minor modifications to wrapper scripts.

## Longer list of differences between SemiBin2 and SemiBin1

The biggest different is that the default training mode is self-supervised mode.

- Output bins are now in a directory called `output_bins` (in
  _SemiBin1_, it actually depended on which parameters were used).
- Output filenames are now anvi'o compatible (effectively, the default value of
  `--tag-output` is `SemiBin`), see discussion at
  [#123](https://github.com/BigDataBiology/SemiBin/issues/123).
- `--compression` defaults to `gz` (instead of `none`)
- ORF finder defaults to the `fast-naive` internal ORF finder
- `--write-pre-reclustering-bins` is `False` by default
- To train in semi-supervised mode, you must use the `train_semi` subcommand
  (and there is no `train` subcommand, you must be specific: `train_semi` or
  `train_self`).

A few arguments that were deprecated before are completely removed:
- `--recluster`: it did nothing already as reclustering is default
- `--mode`: Use `--train-from-many`
- `--training-type`: Use `--semi-supervised` to use semi-supervised learning
  (although that is also deprecated)

