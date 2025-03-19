# Running SemiBin with strobealign-aemb

Strobealign-aemb is a fast abundance estimation method for metagenomic binning available [strobelign](https://github.com/ksahlin/strobealign) version 0.13 (or newer).
For technical reasons, we currently (as of SemiBin 2.1) can only use this mode with a minimum of 5 samples.

## Preparing the data for AEMB

You will need to generate abundance files for each sample in an _all-by-all_ manner.
This means that each sample will be mapped to the contigs of the other samples.

We will illustrate the process with a single sample, but this **must be repeated for all samples**.

1. Split the fasta files using the `split_contigs` subcommand:
```bash
mkdir -p aemb_output/sample1
SemiBin2 split_contigs -i sample1_contigs.fna.gz -o aemb_output/sample1
```

After this step, there will be a file `aemb_output/sample1/split_contigs.fna.gz` that contains both the original contigs as well as split versions (which are required to run SemiBin2).

2. Map reads using [strobealign-aemb](https://github.com/ksahlin/strobealign) to generate the abundance information. Note that version 0.13 (or newer) is required
```bash
strobealign --aemb aemb_output/sample1/split_contigs.fna.gz read1.pair.1.fq.gz read1.pair.2.fq.gz -R 6 -t 8 > sample1_sample1.tsv
strobealign --aemb aemb_output/sample1/split_contigs.fna.gz read2.pair.1.fq.gz read2.pair.2.fq.gz -R 6 -t 8 > sample1_sample2.tsv
strobealign --aemb aemb_output/sample1/split_contigs.fna.gz read3.pair.1.fq.gz read3.pair.2.fq.gz -R 6 -t 8 > sample1_sample3.tsv
strobealign --aemb aemb_output/sample1/split_contigs.fna.gz read4.pair.1.fq.gz read4.pair.2.fq.gz -R 6 -t 8 > sample1_sample4.tsv
strobealign --aemb aemb_output/sample1/split_contigs.fna.gz read5.pair.1.fq.gz read5.pair.2.fq.gz -R 6 -t 8 > sample1_sample5.tsv
```

Each of these commands will generate a file with the abundance information for the sample in the format `sample1_sampleX.tsv`.

## Running SemiBin2

Run SemiBin2 for this sample using the `single_easy_bin` subcommand:

```bash
SemiBin2 single_easy_bin -i contig.fa -a sample1_*.tsv -o aemb_output/sample1
```

⚠️: You should be binning the original contigs, **not** the split contigs.

This will generate the bins in the `aemb_output/sample1` directory.

Note that—from SemiBin2's point-of-view—this is still single sample binning even if the abundance information is generated from multiple samples as the assembly is from a single sample.

## A helper script to run all-by-all abundance estimates

The process above must be repeated for all samples, which can quickly become tedious and error-prone.
Here is a helper script using [Jug](https://jug.readthedocs.io/en/latest/) to automate the process for all samples.
It is most helpful if you are using Jug to parallelize the process, but if you remove the `@TaskGenerator` decorators, it will run sequentially.

It expects the following file structure:

- `samples/` containing the assembled contigs in the format `sample1_assembled.fna.gz`, `sample2_assembled.fna.gz`, ...
- `clean-reads/` containing the reads in the format `sample1.pair.1.fq.gz`, `sample1.pair.2.fq.gz`, `sample2.pair.1.fq.gz`, `sample2.pair.2.fq.gz`, ...
- `aemb_output/` where the output will be written

```python
from jug import TaskGenerator
from jug.utils import jug_execute
import subprocess
from os import makedirs, path
import yaml

samples = [
        'sample0',
        'sample1',
        'sample2',
        'sample3',
        'sample4',
        'sample5',
        'sample6',
        'sample7',
        ]

STROBEALIGN_THREADS = 8
SEMIBIN_THREADS = STROBEALIGN_THREADS

@TaskGenerator
def generate_inputs(s):
    contigs = f'samples/{s}_assembled.fna.gz'
    if not path.exists(contigs):
        raise IOError(f'Expected contig file {f} (for sample {s})')
    out = f'aemb_output/{s}'
    makedirs(out, exist_ok=True)
    subprocess.check_call(
            ['SemiBin2', 'split_contigs',
             '-i', contigs,
             '-o', out])
    return out

@TaskGenerator
def cross_map(ref_out, ref_s, s):
    f1 = f'clean-reads/{s}.pair.1.fq.gz'
    f2 = f'clean-reads/{s}.pair.2.fq.gz'
    if not path.exists(f1):
        raise IOError(f'Expected reads file {f1} (for sample {s})')
    if not path.exists(f2):
        raise IOError(f'Expected reads file {f2} (for sample {s}). Note that {f1} does exist!')
    ofile = f'aemb_output/{ref_s}/mapped_{s}.tsv'
    with open(ofile, 'wb') as out:
        subprocess.check_call(
            ['strobealign',
                '--aemb', f'{ref_out}/split_contigs.fna.gz',
                f1, f2,
                '-t', str(STROBEALIGN_THREADS),
                '-R', '6'],
            stdout=out)
    return ofile


for s in samples:
    out = generate_inputs(s)
    tsv = []
    for s2 in samples:
        tsv.append(cross_map(out, s, s2))

    sb = jug_execute(
        ['SemiBin2', 'single_easy_bin',
            '--threads', str(SEMIBIN_THREADS),
            '-i', f'samples/{s}_assembled.fna.gz',
            '-a'] + tsv + [
                '-o', f'aemb_output/{s}'])

```


