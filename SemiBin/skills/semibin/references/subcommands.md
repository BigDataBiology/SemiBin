# SemiBin2 Subcommand Reference

## All Subcommands Overview

| Subcommand | Purpose | Status |
|---|---|---|
| `single_easy_bin` | Bin contigs (single-sample or co-assembly) in one command | **Primary** |
| `multi_easy_bin` | Bin contigs (multi-sample mode) in one command | **Primary** |
| `concatenate_fasta` | Prepare multi-sample FASTA (rename contigs with sample prefix) | Preprocessing |
| `generate_sequence_features_single` | Generate data.csv/data_split.csv for single/co-assembly | Step-by-step |
| `generate_sequence_features_multi` | Generate data.csv/data_split.csv for multi-sample | Step-by-step |
| `train_self` | Train model with self-supervised learning | Step-by-step |
| `bin` / `bin_short` | Bin short-read contigs using trained model | Step-by-step |
| `bin_long` | Bin long-read contigs using trained model | Step-by-step |
| `split_contigs` | Split contigs for strobealign-aemb pipeline | Preprocessing |
| `check_install` | Verify dependencies are installed | Utility |
| `citation` | Print citation information | Utility |
| `update_model` | Update model file to current format | Utility |
| `generate_cannot_links` / `predict_taxonomy` | Generate cannot-link constraints via mmseqs2 | Deprecated |
| `train_semi` | Train model with semi-supervised learning | Deprecated |
| `download_GTDB` | Download GTDB reference database | Deprecated |

---

## Shared Arguments (available across multiple subcommands)

| Argument | Dest | Type | Default | Subcommands |
|---|---|---|---|---|
| `-i` / `--input-fasta` | `contig_fasta` | str | required | most |
| `-o` / `--output` | `output` | str | required | most |
| `-b` / `--input-bam` | `bams` | list[str] | None | single_easy_bin, multi_easy_bin, generate_sequence_features_* |
| `-a` / `--abundance` | `abundances` | list[str] | None | single_easy_bin, multi_easy_bin, generate_sequence_features_* |
| `-p` / `-t` / `--processes` / `--threads` | `num_process` | int | 0 (all) | most |
| `-m` / `--min-len` | `min_len` | int | None (auto) | most |
| `--ratio` | `ratio` | float | 0.05 | most |
| `--engine` | `engine` | str | `auto` | training/binning cmds |
| `--random-seed` | `random_seed` | int | None | training/binning cmds |
| `--verbose` | `verbose1` | flag | False | most |
| `--quiet` / `-q` | `quiet1` | flag | False | most |
| `--tmpdir` | `tmpdir` | str | None | most |
| `--compression` | `output_compression` | str | `gz` | most |
| `-r` / `--reference-db-data-dir` | `GTDB_reference` | str | None | single_easy_bin, multi_easy_bin, deprecated |
| `--ml-threshold` | `ml_threshold` | int | None | most |

**Auto min-len logic**: If `ratio < 0.05` (proportion of bases in contigs 1000-2500 bp), min-len = 1000; else min-len = 2500. Setting `--min-len` overrides this.

---

## single_easy_bin — Specific Arguments

| Argument | Dest | Type | Default | Notes |
|---|---|---|---|---|
| `--environment` / `--habitat` / `--biome` | `environment` | str | None | Pre-trained model name; only for single BAM |
| `--self-supervised` | `self_supervised` | flag | False | Default training method |
| `--semi-supervised` | `semi_supervised` | flag | False | Deprecated; requires mmseqs2 + GTDB |
| `--sequencing-type` | `sequencing_type` | str | `short_read` | `short_read` or `long_read` |
| `--minfasta-kbs` | `minfasta_kb` | int | 200 | Min output bin size (kbp) |
| `--epochs` / `--epoches` | `epoches` | int | 15 | Training epochs |
| `--batch-size` | `batchsize` | int | 2048 | Training batch size |
| `--orf-finder` | `orf_finder` | str | `fast-naive` | `fast-naive` / `prodigal` / `fraggenescan` |
| `--no-recluster` | `no_recluster` → `recluster` | flag | False | Skip reclustering (good for eukaryotes) |
| `--max-edges` | `max_edges` | int | 200 | Max edges per contig in graph |
| `--max-node` | `max_node` | float | 1.0 | Fraction of contigs to bin (0-1) |
| `--write-pre-reclustering-bins` | `write_pre_reclustering_bins` | bool | False | Save pre-recluster bins |
| `--tag-output` | `output_tag` | str | `SemiBin` | Tag prefix for bin files |
| `--depth-metabat2` | `depth_metabat2` | str | None | Use MetaBAT2 depth file instead of BAM |
| `--cannot-name` | `cannot_name` | str | `cannot` | Cannot-link file name prefix |
| `--taxonomy-annotation-table` | `taxonomy_results_fname` | list[str] | None | Pre-computed mmseqs2 TSV [advanced] |
| `--prodigal-output-faa` | `prodigal_output_faa` | str | None | Bypass ORF calling [deprecated] |

---

## multi_easy_bin — Specific Arguments

All arguments from `single_easy_bin` except `--environment` and `--depth-metabat2`.

| Argument | Dest | Type | Default | Notes |
|---|---|---|---|---|
| `-s` / `--separator` | `separator` | str | `:` | Sample/contig name separator |

**Note**: `--environment` is NOT available for multi_easy_bin (no pre-trained multi-sample models).

---

## concatenate_fasta — Arguments

| Argument | Dest | Type | Default | Notes |
|---|---|---|---|---|
| `-i` / `--input-fasta` | `contig_fasta` | list[str] | required | Multiple FASTA files |
| `-o` / `--output` | `output` | str | required | Output directory |
| `-m` / `--min-len` | `min_len` | int | 0 | Discard shorter contigs |
| `-s` / `--separator` | `separator` | str | `:` | Separator between sample/contig names |
| `--compression` | `output_compression` | str | `gz` | Output compression |

**Output**: `<output>/concatenated.fa.gz` (or `.fa` if `--compression none`)

---

## train_self — Arguments

| Argument | Dest | Type | Default |
|---|---|---|---|
| `-o` / `--output` | `output` | str | required |
| `--data` | `data` | list[str] | required |
| `--data-split` | `data_split` | list[str] | required |
| `--train-from-many` | `train_from_many` | flag | False |
| `--epochs` / `--epoches` | `epoches` | int | 15 |
| `--batch-size` | `batchsize` | int | 2048 |
| `--engine` | `engine` | str | `auto` |
| `-p` / `--processes` | `num_process` | int | 0 |
| `--random-seed` | `random_seed` | int | None |

---

## bin / bin_short — Arguments

| Argument | Dest | Type | Default | Notes |
|---|---|---|---|---|
| `-i` / `--input-fasta` | `contig_fasta` | str | required | |
| `-o` / `--output` | `output` | str | required | |
| `--data` | `data` | str | required | data.csv from generate_sequence_features |
| `--model` | `model_path` | str | None | Path to model.pt |
| `--environment` | `environment` | str | None | Pre-trained model name (alternative to --model) |
| `--minfasta-kbs` | `minfasta_kb` | int | 200 | |
| `--no-recluster` | `no_recluster` | flag | False | |
| `--max-edges` | `max_edges` | int | 200 | |
| `--max-node` | `max_node` | float | 1.0 | |
| `--depth-metabat2` | `depth_metabat2` | str | None | |
| `--orf-finder` | `orf_finder` | str | `fast-naive` | |
| `--engine` | `engine` | str | `auto` | |
| `-p` | `num_process` | int | 0 | |
| `--random-seed` | `random_seed` | int | None | |

---

## bin_long — Arguments

Same as `bin_short` except: no `--no-recluster`, `--max-edges`, `--max-node`.

Default `--minfasta-kbs` is 200 kbp.

---

## Output File Structure

### Single-sample (`single_easy_bin`)

```
output/
├── SemiBinRun.log                    # Full run log (check first on failure)
├── data.csv                          # k-mer + coverage features
├── data_split.csv                    # Features for split contigs (training)
├── model.pt                          # Trained model (if self-supervised)
├── contig_bins.tsv                   # contig -> bin assignment map
├── recluster_bins_info.tsv           # bin summary (short reads, default recluster)
│                                     #   (named bins_info.tsv with --no-recluster or long reads)
└── output_bins/                      # FINAL BINS (default, short AND long reads)
    ├── SemiBin_0.fa.gz
    └── SemiBin_1.fa.gz
```

**Important:** the default final-bins directory is `output_bins/`. The directory
`output_recluster_bins/` is produced ONLY when `--write-pre-reclustering-bins` is
passed; in that case the pre-reclustering bins also go to `output_prerecluster_bins/`.

### Multi-sample (`multi_easy_bin`)

```
multi_output/
├── SemiBinRun.log
├── samples/
│   ├── concatenated.fa.gz            # (if created by concatenate_fasta)
│   ├── S1/
│   │   ├── data.csv
│   │   ├── data_split.csv
│   │   ├── model.pt
│   │   └── output_bins/
│   │       ├── SemiBin_0.fa.gz
│   │       └── (recluster_bins_info.tsv lives in the sample dir, not here)
│   ├── S2/ ...
│   └── S3/ ...
└── bins/                             # All bins, prefixed with sample name
    ├── S1_SemiBin_0.fa.gz
    ├── S2_SemiBin_0.fa.gz
    └── ...
```

---

## Environment → Model Name Mapping

| `--environment` value | Habitat |
|---|---|
| `human_gut` | Human gut metagenomes |
| `dog_gut` | Dog gut |
| `ocean` | Ocean microbiome |
| `soil` | Soil |
| `cat_gut` | Cat gut |
| `human_oral` | Human oral cavity |
| `mouse_gut` | Mouse gut |
| `pig_gut` | Pig gut |
| `built_environment` | Indoor surfaces, built environment |
| `wastewater` | Wastewater treatment |
| `chicken_caecum` | Chicken cecum (contributed post-v1.2) |
| `global` | Combined all-habitat model (use when uncertain) |

---

## generate_sequence_features_single — Arguments

| Argument | Dest | Default | Notes |
|---|---|---|---|
| `-i` / `--input-fasta` | `contig_fasta` | required | |
| `-o` / `--output` | `output` | required | |
| `-b` / `--input-bam` | `bams` | None | One or more BAM/CRAM files |
| `-a` / `--abundance` | `abundances` | None | strobealign-aemb TSV files (≥5 required) |
| `--kmer` | `kmer` | False | Output only k-mer features (no BAM needed) |
| `-m` / `--min-len` | `min_len` | None | |
| `-p` / `--processes` | `num_process` | 0 | |
| `--ml-threshold` | `ml_threshold` | None | |
| `--compression` | `output_compression` | `gz` | |

---

## check_install — Arguments

| Argument | Dest | Default | Notes |
|---|---|---|---|
| `--allow-missing-mmseqs2` | `allow_missing_mmseqs2` | True | mmseqs2 only needed for deprecated semi-supervised mode |

Checks: bedtools, hmmsearch, FragGeneScan, prodigal, mmseqs, torch/CUDA
