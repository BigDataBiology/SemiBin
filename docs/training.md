# Training SemiBin models

SemiBin uses a machine learning approach whereby a model is learned (or _trained_) from data prior to it being used for binning.
Models can be learned from the same data that will be binned.
While this is the simplest approach, it is not the only one as this page explains.

## Single-sample vs multi-sample models

Before we get to different ways of training a model, we note that there are two types of models:

1. Single-sample models
2. Multi-sample models

As the name indicates, single-sample models are used for single-sample binning, while multi-sample models are used for multi-sample binning. 
Multi-sample models are always learned for each individual binning task, but single-sample models can be reused from one sample to the next.
This works best if the samples are somewhat related (_e.g._, from the same habitat), but in that case, it can produce the best results while minizing computational costs.
Thus, we now focus on how single-sample models can be trained.

## Once training is done, single-sample models are interchangeable

**Once a model is trained, it does not matter how it was obtained** and models can be used interchangeably.

## Different ways of training a single-sample model

The simplest approach is to train a model on the sample where it is going to be applied.
However, this can be a computationally-costly endeavour, which is why we explored ways to train and reuse models.

If you have a dataset with multiple similar samples (from the same habitat), you can train on a subset of these and then apply the resulting model to each sample one at a time.
Note that _a single-sample model may be trained from more than one sample_!
This may seem confusing at first, but what makes a model single-sample is that it is going to be applied to each sample independently, not how it was trained.
In fact, we showed that it is beneficial to use several samples to train.

All our prebuilt models are trained from several samples.

## Self-supervised vs. semi-supervised learning

In the [original manuscript describing SemiBin1](https://doi.org/10.1038/s41467-022-29843-y), we presented a semi-supervised approach (and SemiBin is even named after it).
Briefly, contigs are annotated taxonomically (where possible) and then these taxonomic labels are used to learn the model.

Starting with version 1.3, SemiBin also supports _self-supervised_ learning.
With this approach, taxonomic annotation (which was the most time consuming step) is not required.
This approach also appears to not require several samples to learn a good model.
We expect that in version 2, self-supervised learning will become the default mode (any [feedback](mailto:luispedro@big-data-biology.org) or [bug report](https://github.com/BigDataBiology/SemiBin/issues) is welcome).

As of December 2022, a new manuscript describing the self-supervised mode (and other improvements) is under preparation (but we are happy to share a pre-preprint to whomever wants to look at a rough draft and see the data, just [email us](mailto:luispedro@big-data-biology.org)).


