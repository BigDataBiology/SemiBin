# Generating the inputs to SemiBin

## For single-sample mode

Starting with a metagenome, you need to generate a contigs file (`contigs.fna`)
and a sorted BAM file (`output.bam`).

1. Assemble it into a contigs FASTA file. In this case, we are using
   [NGLess](https://ngless.embl.de/) to combine FastQ preprocessing &amp;
   assembly (using
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
write(contigs, ofile='contig.fna')
```

2. Map reads to the FASTA file

### Mapping using NGLess

```python
ngless "1.1"
import "samtools" version "1.0"

input = fastq('sample.fq.gz')
mapped = map(input, fafile='expected.fna')

write(samtools_sort(mapped),
    ofile='output.bam')
```

### Mapping using bowtie2

```bash
bowtie2-build -f contig.fna contig.fna -p 16

bowtie2 -q --fr contig.fna -1 reads_1.fq.gz -2 reads_2.fq.gz -S contig.sam -p 64

samtools view -h -b -S contig.sam -o contig.bam -@ 64

samtools view -b -F 4 contig.bam -o contig.mapped.bam -@ 64

samtools sort -m 1000000000 contig.mapped.bam -o contig.mapped.sorted.bam -@ 64

samtools index contig.mapped.sorted.bam
```

