# Role of gene body methylation in acclimatization and adaptation in a basal metazoan
Groves Dixon, Yi Liao, Line K. Bay and Mikhail V. Matz

datasets and scripts

TagSeq, MBD-seq, and amplicon-bisulfite sequencing datasets are available through NCBI Short Read Archive (Project Accession: SRP049522).

The reads were mapped to the Acropora digitifera genome supplemented with sequences from Symbiodinium: https://www.dropbox.com/s/emn3sxryjnnek91/adig_zoox_combo.tgz?dl=0 

Annotations for this genome version are part of this repository, datasets/adig_zoox_combo.fasta.out.gff.gff3.zip

The file scripts/derivingCounts_bedtools.txt described how counts were extracted from the bam files.

Before playing with R scripts, unzip datasets/mbd_counts.zip and datasets/coords.zip
Start with step1_deseq_GBM.R, step2_deseq_GE.R, and MBDscore_calc.R to generate necessary datasets.
