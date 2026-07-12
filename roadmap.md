# SemiBin2 roadmap/planned work

This document is partly a roadmap and partly a list of issues that we would be happy to see contributions for, but are otherwise likely to work on ourselves as time allows.

**Better integration with AI assistants/agents**. We should make it easy for users who are running an AI agent or assistant to write their pipelines to use SemiBin without mistakes. We took [a few steps](https://semibin.readthedocs.io/en/latest/ai-integration/) in version 2.4.0 with the `install-skills` subcommand and building the [NotebookLM](https://notebooklm.google.com/notebook/e7c2bfae-9756-41d1-99a0-77633c626c58), but there is more to do and we are very open to ideas. Ideally, the AI agent would go beyond being only knowledgeable about SemiBin and would be able to reason about downstream and upstream steps in a pipeline.

**Pre-trained model refresh.** The existing predefined environment models were trained against the [GMGCv1](https://www.nature.com/articles/s41586-021-04233-4) catalog in 2022. The field has evolved since then and it would be good to expand the models.

**A `recluster_only` subcommand** that takes already-clustered bins from another tool (MetaBAT2, VAMB, COMEBin) and reruns SemiBin's marker-gene-based recluster step on them. Cheap to ship if the pieces already exist in `cluster.py:recluster_bins`, and would expose a chunk of SemiBin's quality-improvement value to people who don't want to rebin from scratch.

**Converting to `polars`**. According to several discussions on the [GH issues](https://github.com/BigDataBiology/SemiBin/issues), polars would reduce memory usage. This could be done in a piecemeal fashion, starting with the the bottlenecks and eventually converting all code that uses pandas.

**Enabling model reuse in co-binning.** Right now, the co-binning code path always trains a new model from scratch. In principle, it should be possible to reuse a single pre-trained model for all samples in a co-binning run, but it would imply a lot of internal refactoring (namely in how the `is_combined` flag is currently detected use). This is the underlying reason for [issue #224](https://github.com/BigDataBiology/SemiBin/issues/224) and others.

**Change internal model structure** In particular, using [safetensors](https://github.com/BigDataBiology/SemiBin/issues/206) instead of PyTorch's default format. This could be linked to refreshing the pre-trained models, but it would also be useful on its own.

**Allow users to bypass binning for circular contigs**. This is [issue #124](https://github.com/BigDataBiology/SemiBin/issues/124), one of the oldest open issues.

## More research flavoured ideas

From a more research perspective, there are several ideas to explore within the framework of the self-supervised learning approach that SemiBin uses. This would probably be a good fit for a student project, as it would require a fair amount of experimentation and benchmarking rather than a pull rquest.

**Extend SemiBin to non-bacterial genomes**. SemiBin is currently designed for bacterial genomes, but it would be interesting to extend it to other domains, such as eukaryotic genomes or viral genomes. This would require retraining the models and possibly adjusting the feature extraction methods.

**Explore different features for self-supervised learning**. Currently, SemiBin uses k-mer frequencies, but we could explore other features. In particular, with long-read sequencing, it may be possible to use more complex features on longer contigs.

