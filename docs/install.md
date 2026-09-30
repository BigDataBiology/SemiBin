# Install

SemiBin can run on Python 3.10-3.14.

## Install with pixi

The current recommended way to install SemiBin is to use [pixi](https://pixi.sh/). Pixi will use the packages from conda-forge and bioconda to install SemiBin and its dependencies.

### CPU-only mode

Create a `pixi.toml` file with the following content.

```toml
[workspace]
authors = ["Luis Pedro Coelho <luis@luispedro.org>"]
channels = ["conda-forge", "bioconda"]
name = "semibin_install"
platforms = ["linux-64", "osx-64"]
version = "0.1.0"

[tasks]

[dependencies]
semibin = ">=2.5.0,<3"
```

Then, run `pixi install` in the same directory as the `pixi.toml` file to download and install SemiBin2.


### With GPU support

If you want to use SemiBin with GPU, you need to install Pytorch with GPU support as well. Starting with the example above, you need to add `pytorch-gpu` to the `dependencies` section and `cuda` to the `system-requirements` section.

```toml
[workspace]
authors = ["Luis Pedro Coelho <luis@luispedro.org>"]
channels = ["conda-forge", "bioconda"]
name = "semibin_install"
platforms = ["linux-64"]
version = "0.1.0"

[tasks]

[dependencies]
semibin = ">=2.5.0,<3"
pytorch-gpu = "*"

[system-requirements]
cuda = "12.0"
```

## Install with conda

Pixi is now the recommended way to install SemiBin. However, if you prefer to use conda, you can install SemiBin with it

### Simple Mode

```bash
conda create -n SemiBin
conda activate SemiBin
conda install -c conda-forge -c bioconda semibin
```

### With GPU support

To get GPU support, you need to install Pytorch with GPU support as well, following the instructions on the [Pytorch website](https://pytorch.org/get-started/locally/). As of March 2025, it appears that you need to mix `conda` and `pip` and the following commands should work:

```bash
conda create -n SemiBin python pip
conda activate SemiBin
pip install torch torchvision torchaudio
conda install -c conda-forge -c bioconda semibin
```

This will install SemiBin and all its dependencies. Pixi is—however—the recommended way to install SemiBin with GPU support.

## Install from source

You will need the following dependencies:
- [Bedtools](http://bedtools.readthedocs.org/en/latest/), [Hmmer](http://hmmer.org/), [Samtools](https://github.com/samtools/samtools)

Optional dependencies:
- [Prodigal](https://github.com/hyattpd/Prodigal) (only needed if using `--orf-finder prodigal`; the default `fast-naive` ORF finder does not require it)

You can obtain the required dependencies from conda with the following commands

```bash
conda install -c bioconda bedtools hmmer samtools
```

Then, installing should be a simple matter of running:

```bash
pip install .
```

## Shell tab completion

SemiBin supports tab completion of subcommands and options via [argcomplete](https://kislyuk.github.io/argcomplete/). To enable it in bash or zsh, add the following to your shell startup file (e.g., `~/.bashrc`):

```bash
eval "$(register-python-argcomplete SemiBin2)"
```

(See the [argcomplete documentation](https://kislyuk.github.io/argcomplete/#global-completion) for other shells and for enabling global completion.)

### Static completion scripts with shtab

Alternatively, you can generate a static completion script with
[shtab](https://tqdm.github.io/shtab/). Unlike argcomplete, the generated
script does not call Python each time you press <kbd>Tab</kbd>, so completion
is faster, and it also supports fish, tcsh, and PowerShell. The downside is
that the script needs to be regenerated when you upgrade SemiBin (otherwise,
new subcommands and options will not be completed).

shtab is not a dependency of SemiBin, so install it separately (`pip install
shtab` or `conda install -c conda-forge shtab`). Then, generate the script for
your shell with:

```bash
shtab --shell=bash --prog=SemiBin2 SemiBin.main.semibin_arg_parser > SemiBin2.bash
```

If you have `pixi` installed, you can also generate the script without
installing `shtab` in your current environment. SemiBin itself must be
importable in the temporary environment that `pixi exec` creates, so request it
with `--with`:

```bash
pixi exec -c conda-forge -c bioconda --with semibin shtab -u --shell=bash --prog=SemiBin2 SemiBin.main.semibin_arg_parser > SemiBin2.bash
```

(The `-u` flag makes shtab report an error if `SemiBin.main` cannot be
imported; without it, shtab silently produces an empty script.) If you are
working in a source checkout of SemiBin, running `pixi exec shtab -u ...` from
the top-level directory of the repository also works.

The `--prog=SemiBin2` argument is important: it is the name of the command that
the script will register completions for. Replace `bash` with `zsh`, `fish`,
`tcsh`, or `powershell` as appropriate. Then install the script as follows
(each generated script also has a header with the same instructions):

- **bash**: save it as `~/.local/share/bash_completion/SemiBin2` and add `source ~/.local/share/bash_completion/SemiBin2` to `~/.bashrc`.
- **zsh**: save it as `_SemiBin2` in a directory in your `fpath` (e.g., `~/.local/share/zsh_completion/_SemiBin2`) and add `fpath=(~/.local/share/zsh_completion $fpath)` to `~/.zshrc` (before `compinit`).
- **fish**: save it as `~/.config/fish/completions/SemiBin2.fish`.

For example, for zsh:

```bash
mkdir -p ~/.local/share/zsh_completion
shtab --shell=zsh --prog=SemiBin2 SemiBin.main.semibin_arg_parser > ~/.local/share/zsh_completion/_SemiBin2
```

`SemiBin.main.semibin_arg_parser` is a function that builds and returns the
`argparse.ArgumentParser` for the whole `SemiBin2` command (including all
subcommands), so it can also be used with shtab's Python API
(`shtab.complete(parser, shell="bash")`).

## Alternative ways of running SemiBin

If you use one of these pipelines, we ask that you cite both the pipeline author and the [SemiBin manuscript](https://www.nature.com/articles/s41467-022-29843-y) (as well as any other pipeline-embedded tools which contribute to your results).

- [ATLAS](https://metagenome-atlas.github.io/) is a Snakemake-based pipeline for metagenomics, which includes SemiBin (as well as other binners and tools).
- [Galaxy toolshed](https://toolshed.g2.bx.psu.edu/view/iuc/suite_semibin/) also includes a [SemiBin wrapper](https://toolshed.g2.bx.psu.edu/view/iuc/suite_semibin/)


