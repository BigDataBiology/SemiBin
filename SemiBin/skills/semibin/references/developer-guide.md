# SemiBin2 Developer Guide: Adding New Subcommands

## Module Architecture

```
SemiBin/
├── main.py                   # Entry point: parse_args() + main2() dispatch
├── utils.py                  # validate_normalize_args, fasta_iter, load_fasta,
│                             #   possibly_compressed_write, concatenate_fasta
├── fasta.py                  # fasta_iter (handles .gz/.bz2/.xz transparently)
├── atomicwrite.py            # atomic_write context manager
├── generate_coverage.py      # BAM → per-contig coverage features
├── generate_kmer.py          # k-mer frequency feature computation
├── cluster.py                # Short-read binning/clustering algorithm
├── long_read_cluster.py      # Long-read binning algorithm
├── semi_supervised_model.py  # Semi-supervised model (deprecated)
├── self_supervised_model.py  # Self-supervised model (current default)
├── semibin_version.py        # __version__
└── citation.py               # BIBTEX, RIS, CHICAGO strings
```

**Flow**: `main2(raw_args)` → `parse_args(raw_args)` → `validate_normalize_args()` → dispatch `elif args.cmd == 'X'` → implementation function

---

## The 3-Step Pattern for New Subcommands

### Step 1: Register in `parse_args()` — `SemiBin/main.py`

#### 1a. Create the subparser

Add after the existing `subparsers.add_parser(...)` calls (around line 53–175):

```python
my_cmd = subparsers.add_parser('my_cmd',
                               help='One-line description shown in --help.')
```

#### 1b. Add subcommand-specific arguments

Add unique arguments immediately after the `add_parser` call:

```python
my_cmd.add_argument('-x', '--my-arg',
                    required=True,
                    type=str,
                    help='Description of my argument.',
                    dest='my_arg',        # always use underscores in dest
                    default=None)

my_cmd.add_argument('--my-flag',
                    required=False,
                    help='Enable some feature.',
                    dest='my_flag',
                    action='store_true')
```

#### 1c. Add to shared argument loops

Find the `for p in [...]` loops and add `my_cmd` where appropriate:

```python
# Mandatory -i/--input-fasta and -o/--output (lines ~268-283)
for p in [single_easy_bin, multi_easy_bin, ..., my_cmd]:

# Parallelism: -p/--processes (lines ~381-389)
for p in [train_semi, generate_cannot_links, binning, single_easy_bin,
          multi_easy_bin, ..., training_self, binning_long, my_cmd]:

# --engine for GPU/CPU (lines ~563-569)
for p in [single_easy_bin, multi_easy_bin, train_semi,
          binning, training_self, binning_long, my_cmd]:

# --verbose/--quiet (lines ~438-455)
for p in [single_easy_bin, multi_easy_bin, ..., my_cmd]:

# --random-seed (lines ~546-553)
for p in [train_semi, binning, single_easy_bin, multi_easy_bin,
          training_self, binning_long, my_cmd]:

# --tmpdir (lines ~371-378)
for p in [single_easy_bin, multi_easy_bin, ..., my_cmd]:

# --min-len and --ratio (lines ~399-417)
for p in [single_easy_bin, multi_easy_bin, ..., my_cmd]:

# --compression (lines ~326-338)
for p in [single_easy_bin, ..., concatenate_fasta, my_cmd]:
```

Only add to loops where the argument makes sense for your subcommand.

---

### Step 2: Implement the Function

Implement in `main.py` (for smaller commands) or a new module (for larger ones):

```python
def my_cmd_function(logger, args):
    """
    Implements the my_cmd subcommand.

    Parameters
    ----------
    logger : logging.Logger
        SemiBin2 logger
    args : argparse.Namespace
        Parsed arguments from parse_args()
    """
    import os
    from .atomicwrite import atomic_write
    from .fasta import fasta_iter
    from .utils import possibly_compressed_write

    logger.info('Starting my_cmd...')
    logger.debug(f'Input FASTA: {args.contig_fasta}')

    # Output directory is already created by main2() before dispatch
    output_dir = args.output

    # --- Iterate FASTA file (handles .gz/.bz2/.xz transparently) ---
    results = []
    for header, seq in fasta_iter(args.contig_fasta):
        logger.debug(f'Processing contig: {header} ({len(seq)} bp)')
        results.append((header, len(seq)))

    # --- Write output atomically (safe against partial writes) ---
    output_file = os.path.join(output_dir, 'my_result.tsv')
    with atomic_write(output_file, overwrite=True) as f:
        f.write('contig\tlength\n')
        for header, length in results:
            f.write(f'{header}\t{length}\n')

    # --- Write compressed output ---
    gz_output = os.path.join(output_dir, 'my_output.fa.gz')
    with possibly_compressed_write(gz_output) as f:
        for header, seq in fasta_iter(args.contig_fasta):
            f.write(f'>{header}\n{seq}\n')

    logger.info(f'my_cmd complete. Results in: {output_file}')
```

#### Key conventions:

| Convention | Details |
|---|---|
| Logging | `logger.info()` for user-facing progress; `logger.debug()` for details; `logger.error()` for failures |
| Never use `print()` | Use logger; `print()` is reserved for citation output only |
| File output | Always use `atomic_write()` — prevents corrupt partial files on failure |
| FASTA reading | `fasta_iter(path)` — handles gz/bz2/xz transparently, yields `(header, seq)` tuples |
| Compressed output | `possibly_compressed_write(path)` — writes gz/bz2/xz based on file extension |
| Error exit | `sys.stderr.write(f"Error: ...\n"); sys.exit(1)` OR `logger.error(...); sys.exit(1)` |
| Output dir | `main2()` calls `os.makedirs(args.output, exist_ok=True)` before dispatch — it already exists |

#### Multiprocessing pattern (from `main.py`):

```python
# Pool is defined at module level in main.py:
# Pool = mp.get_context('spawn').Pool

with Pool(min(args.num_process, len(items))) as pool:
    results = [
        pool.apply_async(
            worker_function,
            args=(item, args.output, logger)
        )
        for item in items
    ]
    for r in results:
        result = r.get()  # blocks; raises exception if worker failed
        logger.info(f'Processed: {result}')
```

---

### Step 3: Add Dispatch in `main2()` — `SemiBin/main.py`

In the large `if/elif` dispatch block (starting around line 1523):

```python
elif args.cmd == 'my_cmd':
    my_cmd_function(logger, args)
```

**Important**: Check what `main2()` sets up before dispatch:
- `os.makedirs(args.output, exist_ok=True)` — output dir exists
- `logging.FileHandler` for `SemiBinRun.log` — already set up
- `validate_normalize_args(logger, args)` — args are validated/normalized
- `device` — set for torch subcommands (`args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train_semi', 'bin', 'train_self', 'bin_long']`)
- `contig_dict`, `binned_length` — set for commands in the `generate_cannot_links/generate_sequence_features_single/bin/single_easy_bin/bin_long` group
- CRAM→BAM conversion — done if `hasattr(args, 'bams') and args.bams is not None`

If your command needs GPU, add it to the `args.cmd in [...]` check around line 1477.

---

## Important Argument Naming Conventions

| Argument flag | `dest` | Notes |
|---|---|---|
| `-i`/`--input-fasta` | `contig_fasta` | Used consistently across all subcommands |
| `-b`/`--input-bam` | `bams` | Plural; always a list |
| `-a`/`--abundance` | `abundances` | Plural; always a list |
| `-o`/`--output` | `output` | Output directory path |
| `-p`/`-t`/`--processes`/`--threads` | `num_process` | 0 means all CPUs |
| `--epochs`/`--epoches` | `epoches` | Historical typo kept for backward compat |
| `--batch-size` | `batchsize` | No underscore |
| `--min-len` | `min_len` | |
| `--random-seed` | `random_seed` | |
| `--no-recluster` | `no_recluster` | In `parse_args()`: inverted to `args.recluster = not args.no_recluster` |
| `--sequencing-type` | `sequencing_type` | |
| `--minfasta-kbs` | `minfasta_kb` | Note: kbs → kb (singular) |
| `--reference-db-data-dir` | `GTDB_reference` | All caps |
| `-s`/`--separator` | `separator` | Multi-sample only |

---

## `validate_normalize_args()` Hook — `SemiBin/utils.py`

This function runs before dispatch and handles:
- Checking that input files exist (using `expect_file()` / `expect_file_list()`)
- Setting `args.num_process` to CPU count if 0
- Calling `check_training_type()` to set `args.training_type`
- Validating `args.engine`

If your subcommand has inputs that need validation, add them here:

```python
# In validate_normalize_args(), find the relevant hasattr(args, ...) blocks:
if hasattr(args, 'my_arg'):
    expect_file(args.my_arg)
```

---

## Verbose/Quiet Argument Hack

The `--verbose`/`--quiet` arguments use a special workaround because argparse can't accept the same option in both global and subcommand scope cleanly.

- Subcommand parsers use `dest='verbose1'` and `dest='quiet1'`
- At the end of `parse_args()`, a loop maps `verbose1`→`verbose` and `quiet1`→`quiet`
- When you add your parser to the `--verbose/--quiet` loop, use the standard `dest='verbose1'` / `dest='quiet1'` (this is already handled by the loop you add to)

---

## GPU/Device Selection

Commands that use PyTorch training are listed in this check in `main2()` (around line 1477):

```python
if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train_semi', 'bin', 'train_self', 'bin_long']:
    import torch
    if args.engine == 'cpu':
        device = torch.device("cpu")
    else:
        if torch.cuda.is_available():
            device = torch.device("cuda")
        else:
            device = torch.device("cpu")
            logger.warning('Did not detect GPU ...')
```

Add your command to this list if it uses PyTorch. The `device` variable is then available in the dispatch block.

---

## Alias Pattern

To create an alias for your command (like `bin` / `bin_short`):

```python
# In parse_args():
my_cmd = subparsers.add_parser('my_cmd', aliases=['my_cmd_alias'])

# In main2(), normalize the alias:
if args.cmd == 'my_cmd_alias':
    args.cmd = 'my_cmd'
```

---

## Testing

```python
# Programmatic testing via main2():
from SemiBin.main import main2

# Basic invocation test:
main2(['my_cmd', '-i', 'test_contigs.fa', '-o', 'test_output'])

# With optional args:
main2(['my_cmd', '-i', 'contigs.fa', '-o', 'output',
       '--my-arg', 'value', '--verbose'])

# Test error handling:
import pytest
with pytest.raises(SystemExit):
    main2(['my_cmd'])  # Missing required args
```

---

## Common Patterns from Existing Code

### Loading and filtering a FASTA file

```python
from .utils import load_fasta
# Returns: (c_min_len, must_link_threshold, contig_dict)
# contig_dict: {header: sequence_string}
c_min_len, must_link_threshold, contig_dict = load_fasta(args.contig_fasta, args.ratio)
binned_length = c_min_len if args.min_len is None else args.min_len
```

### Checking for empty or insufficient input

```python
if not contig_dict:
    logger.error(f'Input file {args.contig_fasta} is empty.')
    sys.exit(1)

n_pass = sum(len(c) >= binned_length for c in contig_dict.values())
if n_pass < 4:
    logger.error(f'Only {n_pass} contig(s) above minimum length {binned_length} bp.')
    sys.exit(1)
```

### Subprocess execution

```python
import subprocess
try:
    subprocess.check_call(['tool', 'arg1', 'arg2'], stdout=None)
except Exception as e:
    sys.stderr.write(f"Error: Running tool failed (error: {e})\n")
    sys.exit(1)
```

### Writing gzip-compressed FASTA

```python
from .utils import possibly_compressed_write

output_fa = os.path.join(args.output, 'output.fa.gz')
with possibly_compressed_write(output_fa) as out:
    for header, seq in fasta_iter(args.contig_fasta):
        out.write(f'>{header}\n{seq}\n')
```
