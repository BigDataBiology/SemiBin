# AGENTS.md

## Project Overview

SemiBin is a metagenomic binning tool that uses siamese neural networks (PyTorch) to cluster contigs from metagenomic assemblies into genome bins. It supports both self-supervised (SemiBin2, default) and semi-supervised learning modes, with pre-trained models for 11 environments.

## Common Commands

### Running tests
```bash
python -m pytest                          # all unit tests
python -m pytest test/test_bin.py         # single test file
python -m pytest test/test_bin.py -k test_name  # single test
```

### Integration tests (require installed package + external tools)
```bash
uv pip install .
python integration-tests/easy_commands2.py   # run one integration test
```

### Environment setup (pixi manages conda deps including external tools)
```bash
pixi run -e test-py312 pytest    # run tests in a specific Python version env
```

### Install for development
```bash
uv pip install -e .
```

## Architecture

### Entry point and CLI
`SemiBin/main.py` — CLI argument parsing (`parse_args`) and all subcommand handlers. The entry point is `SemiBin.main:main2`, exposed as `SemiBin2` command. This is a large file (~1600 lines) that orchestrates the entire pipeline.

### Neural network models
- `SemiBin/semi_supervised_model.py` — Siamese network architecture (`Semi_encoding_single`, `Semi_encoding_multiple`) and training loop using must-link/cannot-link constraints
- `SemiBin/self_supervised_model.py` — Contrastive learning training (`train_self`), reuses the network classes from `semi_supervised_model.py`

### Feature generation
- `SemiBin/generate_coverage.py` — Converts BAM alignments to coverage/abundance features via samtools+bedtools
- `SemiBin/generate_kmer.py` — Extracts k-mer composition features from FASTA sequences
- `SemiBin/markers.py` — HMM marker gene detection using HMMER (bundled `marker.hmm` profile)

### Clustering and binning
- `SemiBin/cluster.py` — Main clustering via `run_embed_infomap()` (graph-based with python-igraph) and `recluster_bins()` for quality refinement
- `SemiBin/long_read_cluster.py` — DBSCAN-based ensemble clustering for long reads

### Supporting modules
- `SemiBin/utils.py` — Shared utilities: data normalization, constraint generation, FASTA loading, bin writing
- `SemiBin/atomicwrite.py` — Crash-safe file writing
- `SemiBin/orffinding.py` / `naive_orffinder.py` — ORF prediction (prodigal/fraggenescan/built-in)
- `SemiBin/fasta.py` — FASTA parsing
- `SemiBin/gtdb.py` — GTDB reference database integration

### Pre-trained models
`SemiBin/models/*.pt` — 11 environment-specific models (human_gut, ocean, soil, etc.)

### Key data flow
1. **Feature generation**: contigs.fa + BAM → coverage features + k-mer features
2. **Training**: features → self-supervised contrastive learning (or semi-supervised with taxonomy constraints)
3. **Clustering**: trained model embeds contigs → infomap graph clustering (or DBSCAN for long reads)
4. **Output**: clustered bins written as FASTA files

## External tool dependencies
Managed via conda/pixi (not pip): samtools, bedtools, hmmer, mmseqs2, prodigal, fraggenescan.

## Commit message conventions
Prefix commits with: `BUG`, `ENH`, `MIN`, `RFCT`, `TST`, `DOC`, `RLS`. Can be combined (e.g., `BUG+DOC`).

## Changelog and release notes
For new features and bugfixes, update both:
- `ChangeLog` — add a line under `Unreleased` summarising the change (include issue number if applicable)
- `docs/whatsnew.md` — add a bullet under the current unreleased version section

## Key details
- Python 3.8+ compatibility required
- Uses `mp.get_context('spawn').Pool` for multiprocessing (not fork)
- `hypothesis` pinned to `<=6.112.1` in test dependencies
