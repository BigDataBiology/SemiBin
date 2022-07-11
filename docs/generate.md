# Generating the inputs to SemiBin (from a metagenome)

Starting with a metagenome, you need to generate a contigs file (`contigs.fa`)
and a sorted BAM file (`output.bam`) from mapping the metagenomic reads to the
assembled contigs.

**Step 1**: Assemble it into a contigs FASTA file. In this case, we are using
[NGLess](https://ngless.embl.de/) to combine FastQ preprocessing &amp; assembly
(using
[MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884),
but any other system will work.

```python
ngless "1.2"

input = paired('reads_1.fq.gz', 'reads_2.fq.gz')
input = preprocess(input) using |r|:
    r = substrim(r, min_quality=25)
    if len(r) < 45:
        discard

contigs = assemble(input)
write(contigs, ofile='contig.fa')
```

**Step 2**: Map reads to the FASTA file generated in _Step 1_.

### Mapping using NGLess

```python
ngless "1.2"
import "samtools" version "1.0"

input = fastq('sample.fq.gz')
mapped = map(input, fafile='expected.fa')

write(samtools_sort(mapped),
    ofile='output.bam')
```

### Mapping using bowtie2

You can also use `bowtie2` directly, for example (using 4 threads, you can
adjust `-p 4` as needed when calling `bowtie2`):

```bash
bowtie2-build -f contig.fa contig.fa

bowtie2 -q --fr -x contig.fa -1 reads_1.fq.gz -2 reads_2.fq.gz -S contig.sam -p 4

samtools view -h -b -S contig.sam -o contig.bam
samtools view -b -F 4 contig.bam -o contig.mapped.bam
samtools sort contig.mapped.bam -o contig.mapped.sorted.bam

samtools index contig.mapped.sorted.bam
```

### Generate cannot-link constraints using CAT (advanced)

_**Note**: Unless you understand exactly what is going on, you probably **do not** want to do this. Feel free to check in [with us](https://groups.google.com/g/semibin-users) if you have doubts._

SemiBin uses mmseqs2 by default, but you can also use [CAT](https://github.com/dutilh/CAT) to produce contig taxonomic classifications and generate the cannot-link pairs.

```bash
CAT contigs \
        -c contig.fa \
        -d CAT_prepare_20200304/2020-03-04_CAT_database \
        --path_to_prodigal $path_to_prodigal \
        --path_to_diamond $path_to_diamond \
        -t CAT_prepare_20200304/2020-03-04_taxonomy \
        -o CAT_output/CAT \
        --force \
        -f 0.5 \
        --top 11 \
        --I_know_what_Im_doing \
        --index_chunks 1

CAT add_names \
    CAT_output/CAT.contig2classification.txt \
    -o CAT_output/CAT.out \
    -t CAT_prepare_20200304/2020-03-04_taxonomy \
    --force \
    --only_official
```

Generate cannot-link constrains using the script `generate_cannot_link.py` in the `scripts/` directory

```bash
python script/generate_cannot_link.py -i CAT.out -c contig.fa -s sample-name -o output --CAT
```

