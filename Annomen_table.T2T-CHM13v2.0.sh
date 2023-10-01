#!/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

# Requirements
# 1. Perl: https://www.perl.org
# 2. BioPerl: http://www.bioperl.org/wiki/Main_Page
#    - Bio::DB::Fasta
#    - Bio::SeqIO
# 3. EMBOSS: http://emboss.sourceforge.net or EMBOSS-6.6.0.reduced.tar.gz or conda install -c bioconda emboss
#    - needle
#    - stretcher
# 4. Basic linux commands: bash, rm, gzip, sort, echo, find, sed, awk, wget

# Remove old files
rm -rf refseq refseq.transcript.fasta refseq.protein.fasta refseq.transcript.gbff
rm -rf Annomen_table.txt Annomen_table.log

# Download RefSeq files from NCBI FTP
wget --no-verbose --no-check-certificate --recursive --execute robots=off --reject 'index.html*' --no-parent --level=1 --no-host-directories --cut-dirs=0 https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz;     do gzip -dc $file; done > refseq.transcript.fasta
for file in refseq/H_sapiens/mRNA_Prot/human.*.protein.faa.gz; do gzip -dc $file; done > refseq.protein.fasta
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - refseq.transcript.gbff

# Download genome files from NCBI FTP
wget --no-verbose --no-check-certificate --recursive --execute robots=off --reject 'index.html*' --no-parent --level=1 --no-host-directories --cut-dirs=6 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/

gzip -dc GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz | sed 's/^>.* chromosome />chr/' | sed 's/, .*$//' > genome.fasta

grep -v '^#' GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt | cut -f7,10 | sed 's/\r$//' > chromosome.UCSC.txt

# Generate Annomen table
perl Annomen_table.pl -c chromosome.UCSC.txt GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz genome.fasta refseq.transcript.fasta refseq.protein.fasta refseq.transcript.gbff > Annomen_table.txt 2> Annomen_table.log
