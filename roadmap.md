# SemiBin2 roadmap/planned work

This document is partly a roadmap and partly a list of issues that we would be happy to see contributions for.

**Better integration with AI assistants/agents**. We should make it easy for users who are running an AI agent or assistant to write their pipelines to use SemiBin without mistakes.

**Pre-trained model refresh.** The existing predefined environment models were trained against the [GMGC](https://www.nature.com/articles/s41586-021-04233-4) catalog in 2022. The field has evolved since then and it would be good to expand the models.

**A `recluster_only` subcommand** that takes already-clustered bins from another tool (MetaBAT2, VAMB, COMEBin) and reruns SemiBin's marker-gene-based recluster step on them. Cheap to ship if the pieces already exist in `cluster.py:recluster_bins`, and would expose a chunk of SemiBin's quality-improvement value to people who don't want to rebin from scratch.

**Converting to `polars`**. According to several discussions on the [GH issues](https://github.com/BigDataBiology/SemiBin/issues), polars would reduce memory usage.

**Enabling model reuse in co-binning.** Right now, the co-binning code path always trains a new model from scratch. In principle, it should be possible to reuse a single pre-trained model for all samples in a co-binning run, but it would imply a lot of internal refactoring (namely in how the `is_combined` flag is currently detected use).

**Change internal model structure** In particular, using [safetensors](https://github.com/BigDataBiology/SemiBin/issues/206) instead of PyTorch's default format. This could be linked to refreshing the pre-trained models, but it would also be useful on its own.

**Allow users to bypass binning for circular contigs**. This is [issue #124](https://github.com/BigDataBiology/SemiBin/issues/124), one of the oldest open issues.


