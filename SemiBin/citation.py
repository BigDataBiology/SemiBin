BIBTEX = '''\
@ARTICLE{Pan2022semibin,
  title     = "A deep siamese neural network improves metagenome-assembled
               genomes in microbiome datasets across different environments",
  author    = "Pan, Shaojun and Zhu, Chengkai and Zhao, Xing-Ming and Coelho,
               Luis Pedro",
  abstract  = "Metagenomic binning is the step in building metagenome-assembled
               genomes (MAGs) when sequences predicted to originate from the
               same genome are automatically grouped together. The most
               widely-used methods for binning are reference-independent,
               operating de novo and enable the recovery of genomes from
               previously unsampled clades. However, they do not leverage the
               knowledge in existing databases. Here, we introduce SemiBin, an
               open source tool that uses deep siamese neural networks to
               implement a semi-supervised approach, i.e. SemiBin exploits the
               information in reference genomes, while retaining the capability
               of reconstructing high-quality bins that are outside the
               reference dataset. Using simulated and real microbiome datasets
               from several different habitats from GMGCv1 (Global Microbial
               Gene Catalog), including the human gut, non-human guts, and
               environmental habitats (ocean and soil), we show that SemiBin
               outperforms existing state-of-the-art binning methods. In
               particular, compared to other methods, SemiBin returns more
               high-quality bins with larger taxonomic diversity, including
               more distinct genera and species.",
  journal   = "Nat. Commun.",
  publisher = "Nature Publishing Group",
  volume    =  13,
  number    =  1,
  pages     = "2326",
  month     =  apr,
  year      =  2022,
  language  = "en",
  doi       = "10.1038/s41467-022-29843-y"
}

@ARTICLE{Pan2023semibin2,
  title    = "{SemiBin2}: self-supervised contrastive learning leads to better
              {MAGs} for short- and long-read sequencing",
  author   = "Pan, Shaojun and Zhao, Xing-Ming and Coelho, Luis Pedro",
  abstract = "MOTIVATION: Metagenomic binning methods to reconstruct
              metagenome-assembled genomes (MAGs) from environmental samples
              have been widely used in large-scale metagenomic studies. The
              recently proposed semi-supervised binning method, SemiBin,
              achieved state-of-the-art binning results in several
              environments. However, this required annotating contigs, a
              computationally costly and potentially biased process. RESULTS:
              We propose SemiBin2, which uses self-supervised learning to learn
              feature embeddings from the contigs. In simulated and real
              datasets, we show that self-supervised learning achieves better
              results than the semi-supervised learning used in SemiBin1 and
              that SemiBin2 outperforms other state-of-the-art binners.
              Compared to SemiBin1, SemiBin2 can reconstruct 8.3-21.5\% more
              high-quality bins and requires only 25\% of the running time and
              11\% of peak memory usage in real short-read sequencing samples.
              To extend SemiBin2 to long-read data, we also propose
              ensemble-based DBSCAN clustering algorithm, resulting in
              13.1-26.3\% more high-quality genomes than the second best binner
              for long-read data. AVAILABILITY AND IMPLEMENTATION: SemiBin2 is
              available as open source software at
              https://github.com/BigDataBiology/SemiBin/ and the analysis
              scripts used in the study can be found at
              https://github.com/BigDataBiology/SemiBin2\_benchmark.",
  journal  = "Bioinformatics",
  volume   =  39,
  number   = "39 Suppl 1",
  pages    = "i21--i29",
  month    =  jun,
  year     =  2023,
  language = "en",
  doi      = "10.1093/bioinformatics/btad209"
}
'''

RIS = '''\
TY  - JOUR
AU  - Pan, Shaojun
AU  - Zhu, Chengkai
AU  - Zhao, Xing-Ming
AU  - Coelho, Luis Pedro
PY  - 2022
DA  - 2022/04/28
TI  - A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments
JO  - Nature Communications
SP  - 2326
VL  - 13
IS  - 1
AB  - Metagenomic binning is the step in building metagenome-assembled genomes (MAGs) when sequences predicted to originate from the same genome are automatically grouped together. The most widely-used methods for binning are reference-independent, operating de novo and enable the recovery of genomes from previously unsampled clades. However, they do not leverage the knowledge in existing databases. Here, we introduce SemiBin, an open source tool that uses deep siamese neural networks to implement a semi-supervised approach, i.e. SemiBin exploits the information in reference genomes, while retaining the capability of reconstructing high-quality bins that are outside the reference dataset. Using simulated and real microbiome datasets from several different habitats from GMGCv1 (Global Microbial Gene Catalog), including the human gut, non-human guts, and environmental habitats (ocean and soil), we show that SemiBin outperforms existing state-of-the-art binning methods. In particular, compared to other methods, SemiBin returns more high-quality bins with larger taxonomic diversity, including more distinct genera and species.
SN  - 2041-1723
UR  - https://doi.org/10.1038/s41467-022-29843-y
DO  - 10.1038/s41467-022-29843-y
ID  - Pan2022
ER  - 

TY  - JOUR
AU  - Pan, Shaojun
AU  - Zhao, Xing-Ming
AU  - Coelho, Luis Pedro
T1  - SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing
PY  - 2023
Y1  - 2023/06/01
DO  - 10.1093/bioinformatics/btad209
JO  - Bioinformatics
JA  - Bioinformatics
VL  - 39
IS  - Supplement_1
SP  - i21
EP  - i29
SN  - 1367-4811
AB  - Metagenomic binning methods to reconstruct metagenome-assembled genomes (MAGs) from environmental samples have been widely used in large-scale metagenomic studies. The recently proposed semi-supervised binning method, SemiBin, achieved state-of-the-art binning results in several environments. However, this required annotating contigs, a computationally costly and potentially biased process.We propose SemiBin2, which uses self-supervised learning to learn feature embeddings from the contigs. In simulated and real datasets, we show that self-supervised learning achieves better results than the semi-supervised learning used in SemiBin1 and that SemiBin2 outperforms other state-of-the-art binners. Compared to SemiBin1, SemiBin2 can reconstruct 8.3–21.5% more high-quality bins and requires only 25% of the running time and 11% of peak memory usage in real short-read sequencing samples. To extend SemiBin2 to long-read data, we also propose ensemble-based DBSCAN clustering algorithm, resulting in 13.1–26.3% more high-quality genomes than the second best binner for long-read data.SemiBin2 is available as open source software at https://github.com/BigDataBiology/SemiBin/ and the analysis scripts used in the study can be found at https://github.com/BigDataBiology/SemiBin2_benchmark.
Y2  - 10/22/2023
UR  - https://doi.org/10.1093/bioinformatics/btad209
ER  - 
'''


CHICAGO = '''\
Pan, Shaojun, Chengkai Zhu, Xing-Ming Zhao, and Luis Pedro Coelho. 2022. "A Deep Siamese Neural Network Improves Metagenome-Assembled Genomes in Microbiome Datasets across Different Environments." Nature Communications 13 (1): 2326. https://doi.org/10.1093/bioinformatics/btad209

Pan, Shaojun, Xing-Ming Zhao, and Luis Pedro Coelho. 2023. "SemiBin2: Self-Supervised Contrastive Learning Leads to Better MAGs for Short- and Long-Read Sequencing." Bioinformatics  39 (39 Suppl 1): i21–29. https://doi.org/10.1038/s41467-022-29843-y
'''
