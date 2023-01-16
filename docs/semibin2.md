# SemiBin2

Starting with version 1.5 (officially _SemiBin2 beta_), installing the SemiBin
package installs two scripts: `SemiBin` and `SemiBin2`.

They have the same functionality, but slightly different interfaces. The exact
interface to `SemiBin2` should be considered as unstable (while we will strive
to maintain backwards compatibility if you call the `SemiBin` script).

# Differences between SemiBin2 and SemiBin1

The biggest different is that the default training mode is self-supervised mode.

- Output bins are now **always** in a directory called `output_bins` (in
- Output filenames are now anvi'o compatible (effectively, the default value of `--tag-output` is `SemiBin`) (see discussion in [#123](https://github.com/BigDataBiology/SemiBin/issues/123))
  _SemiBin1_, it actually depended on which parameters were used)
- `--compression` defaults to `gz` (instead of `none`)
- ORF finder defaults to the `fast-naive` internal ORF finder
- `--write-pre-reclustering-bins` is `False` by default
- To train in semi-supervised mode, you must use the `train_semi` subcommand
  (and there is no `train` subcommand, you must be specific: `train_semi` or
  `train_self`).

A few arguments that were deprecated before are completely removed:
- `--recluster`: it did nothing already as reclustering is default
- `--mode`: Use `--train-from-many`
- `--training-type`: Use `--semi-supervised` to use semi-supervised learning (although that is also deprecated)

