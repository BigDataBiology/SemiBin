---
name: semibin
description: >
  Auto-activate when the user mentions SemiBin, SemiBin2, metagenomic binning, MAGs
  (metagenome-assembled genomes), binning contigs, or asks about tools like MetaBAT2,
  MaxBin, CONCOCT, or DAS Tool in the context of comparing to or replacing them with
  SemiBin2. Also activate when the user is working with BAM files alongside metagenomic
  assemblies, or discusses single-sample/multi-sample/co-assembly binning workflows.
version: 1.0.0
---

# Using SemiBin2 Correctly

SemiBin2 is a neural-network metagenomic binner: it groups assembled contigs into
bins (putative genomes / MAGs). The command is `SemiBin2`. This skill exists so an
agent can construct a *correct* invocation on the first try and read the output from
the right place.

## Decision tree: which command do I run?

Ask two questions — (a) how many assemblies, and (b) do you have a matching
pre-trained model — then pick:

- **One assembly, one sample, habitat has a built-in model** → `single_easy_bin --environment <habitat>` (fastest; minutes).
- **One assembly, one sample, no matching model** → `single_easy_bin --self-supervised`.
- **One co-assembly, reads from several samples** → `single_easy_bin --self-supervised` with *multiple* `-b` BAMs (pre-trained models CANNOT be used with >1 BAM).
- **Several assemblies binned together (multi-sample)** → `concatenate_fasta` first, then `multi_easy_bin` (no pre-trained models).
- **Long-read / hybrid assembly** → add `--sequencing-type long_read` to any of the above (hybrid is treated as long_read).

`single_easy_bin` and `multi_easy_bin` are the one-shot commands and are almost
always what you want. The individual steps (`generate_sequence_features_*`,
`train_self`, `bin_short`/`bin_long`) exist only for splitting work across a cluster;
prefer the easy commands unless the user explicitly wants staged execution.

## Copy-paste command templates

Single-sample, known habitat (fastest path):
```bash
SemiBin2 single_easy_bin --environment human_gut -i contigs.fa -b sample.sorted.bam -o output
```

Single-sample, unknown habitat (train a new model):
```bash
SemiBin2 single_easy_bin --self-supervised -i contigs.fa -b sample.sorted.bam -o output
```

Co-assembly (one FASTA, several BAMs):
```bash
SemiBin2 single_easy_bin --self-supervised -i coassembly.fa -b S1.sorted.bam S2.sorted.bam S3.sorted.bam -o output
```

Multi-sample (two steps):
```bash
SemiBin2 concatenate_fasta -i S1.fa S2.fa S3.fa -o concat_out
# Map each sample's reads to concat_out/concatenated.fa.gz, producing sorted+indexed BAMs
SemiBin2 multi_easy_bin -i concat_out/concatenated.fa.gz -b S1.sorted.bam S2.sorted.bam S3.sorted.bam -o multi_out
```

Long reads: append `--sequencing-type long_read` to any command above.

## Inputs (validate these before running)

- **FASTA of contigs** — gzip/bzip2/xz accepted; contig names must be unique; need
  several contigs above the min length (default min-len is auto-chosen: 1000 or 2500 bp).
- **BAM/CRAM** — reads mapped to *that same* FASTA, **sorted** (`samtools sort`) and
  **indexed** (`samtools index`, a `.bai`/`.crai` must exist next to it). Passing `-b`
  multiple files = co-assembly/multi-sample coabundance.
- **Multi-sample concatenated FASTA** — contig headers must be `<sample><sep><contig>`
  (default separator `:`). Do NOT hand-build this; use `SemiBin2 concatenate_fasta`.
- **Alternative to BAMs**: `-a/--abundance` accepts strobealign-aemb abundance TSVs,
  but ONLY with **≥5 samples**.

## Output layout (where the bins actually are)

By default reclustering is ON for short reads. The final bins land in **`output/output_bins/`**:

| Situation | Final bins directory | Summary table |
|---|---|---|
| short reads, default | `output/output_bins/` | `recluster_bins_info.tsv` |
| short reads, `--no-recluster` | `output/output_bins/` | `bins_info.tsv` |
| long reads | `output/output_bins/` | `bins_info.tsv` |
| short reads, `--write-pre-reclustering-bins` | `output/output_recluster_bins/` (final) + `output/output_prerecluster_bins/` | `recluster_bins_info.tsv` |
| multi-sample | `multi_out/bins/` (all, sample-prefixed) and `multi_out/samples/<S>/output_bins/` | per-sample as above |

> The default final directory is `output_bins/`. `output_recluster_bins/` exists ONLY
> when `--write-pre-reclustering-bins` is passed — do not point users there by default.

Other useful files in `output/`:
- `contig_bins.tsv` — contig → bin assignment map.
- `*_bins_info.tsv` — one row per bin: name, total_bp, n_contigs, N50, L50.
- `model.pt` — the trained model (self-supervised runs).
- `SemiBinRun.log` — timestamped run log, including any fatal error. **First place to look when a run fails.**

Bin FASTA files are named `<tag>_<n>.fa.gz` (tag defaults to `SemiBin`, set via `--tag-output`).

## Pre-trained environments

`human_gut` · `dog_gut` · `ocean` · `soil` · `cat_gut` · `human_oral` · `mouse_gut` ·
`pig_gut` · `built_environment` · `wastewater` · `chicken_caecum` · `global`

Use `--environment global` when the habitat is unknown but you still want the fast path.
Built-in models work **only** in single-sample mode (exactly one BAM).

## Flags worth knowing

- `--self-supervised` / `--semi-supervised` — training mode. Self-supervised is the modern default; semi-supervised is deprecated (needs mmseqs2 + GTDB).
- `--sequencing-type short_read|long_read|hybrid` — only on the easy commands; hybrid ⇒ long_read.
- `--no-recluster` — skip the reclustering refinement (recommended for eukaryotic metagenomes; the recluster step is tuned for prokaryotes).
- `--orf-finder fast-naive|prodigal|fraggenescan` — **default is `fast-naive`, which needs no external tool.** Only switch if the user specifically wants prodigal/FragGeneScan.
- `-p/-t/--processes/--threads` — CPUs; `0` (default) means use all.
- `--engine auto|gpu|cpu` — device selection; `auto` uses GPU if present.
- `--random-seed <int>` — set for reproducible runs (default is nondeterministic).
- `--environment` also accepts aliases `--habitat` / `--biome`.

## Common failure modes → fix

1. **"provided pretrained model only used in single-sample binning"** — you passed
   `--environment` with multiple BAMs. Fix: use `--self-supervised` for co-assembly/multi-sample.
2. **"Expected contigs to contain separator character"** — multi-sample FASTA lacks the
   `<sample>:<contig>` prefix. Fix: build it with `SemiBin2 concatenate_fasta`.
3. **Poor bins on long-read/hybrid data** — forgot `--sequencing-type long_read`. Always set it for ONT/PacBio/hybrid.
4. **"abundances from strobealign-aemb can only be used with at least 5 samples"** — `-a` needs ≥5 samples; use `-b` BAMs instead for fewer.
5. **BAM errors / no coverage** — BAM not sorted or not indexed, or mapped to a different FASTA. Re-`samtools sort` + `samtools index` against the exact input FASTA.
6. **Eukaryotic genomes fragmented** — add `--no-recluster`.
7. **"all are shorter than N basepairs" / "only N contig(s) contain at least N basepairs"** — too few contigs pass the length filter (need ≥4 above min-len). Lower `-m/--min-len` (e.g. `500`) or improve the assembly.
8. **GPU/CUDA errors** — force CPU with `--engine cpu`. For "model format" errors on an old model, run `SemiBin2 update_model -m old.pt -o new.pt`.
9. **"Running mmseqs createdb failed"** — only affects the deprecated semi-supervised mode; switch to `--self-supervised` (no mmseqs2/GTDB needed).
10. **Run aborted, unclear why** — read `output/SemiBinRun.log`.

## Sanity checks before you claim success

- Run `SemiBin2 check_install --verbose` if unsure the install is complete.
- After a run, confirm `output/output_bins/` is non-empty and inspect `*_bins_info.tsv`
  for bin count and sizes. A run that logs "No bins were created" needs more/longer
  contigs or better coverage, not a rerun with the same inputs.
- Downstream quality is assessed with **CheckM2** (recommended), CheckM, or GTDB-Tk —
  SemiBin2 does not compute completeness/contamination itself.

## Installation

```bash
pixi global install semibin                              # recommended
conda install -c conda-forge -c bioconda semibin         # conda
SemiBin2 check_install --verbose                          # verify
```

## Generating inputs (assembly + mapping)

If the user has only reads, they must assemble and map first:
```bash
# assemble reads → contig.fa (MEGAHIT/metaSPAdes/etc.), then map:
bowtie2-build -f contig.fa contig.fa
bowtie2 -q --fr -x contig.fa -1 reads_1.fq.gz -2 reads_2.fq.gz -p 4 -S contig.sam
samtools view -bS contig.sam | samtools view -b -F 4 - -o contig.mapped.bam
samtools sort contig.mapped.bam -o sample.sorted.bam
samtools index sample.sorted.bam
```

## References

- `references/subcommands.md` — full per-subcommand argument tables and output trees.
- `references/developer-guide.md` — module map and the pattern for adding a new subcommand (for work *inside* the SemiBin codebase).
