#!/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

# Requirements
# 1. Perl: https://www.perl.org
# 2. BioPerl: http://www.bioperl.org/wiki/Main_Page
#    - Bio::DB::Fasta
#    - Bio::SeqIO
# 3. EMBOSS: http://emboss.sourceforge.net
#    - needle
#    - stretcher
# 4. Basic linux commands: bash, rm, gzip, sort, echo, find, sed, awk, wget

# Remove old files
rm -rf mm10.fa.gz mm10.fa
rm -rf refseq refseq.transcript.fasta refseq.protein.fasta refseq.transcript.gbff
rm -rf Annomen_table.txt Annomen_table.log

# Prepare reference genome fasta file
if test -r "$genomeFastaFile"; then
	ln -sf "$genomeFastaFile" genome.fasta
else
	wget --no-verbose --no-check-certificate https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
	gzip -d mm10.fa.gz
	ln -sf mm10.fa genome.fasta
fi

# Download RefSeq files from NCBI FTP
wget --no-verbose --no-check-certificate --recursive -e robots=off --reject "index.html" --no-host-directories --cut-dirs=0 https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.fna.gz;     do gzip -dc $file; done > refseq.transcript.fasta
for file in refseq/M_musculus/mRNA_Prot/mouse.*.protein.faa.gz; do gzip -dc $file; done > refseq.protein.fasta
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - refseq.transcript.gbff

wget --no-verbose --no-check-certificate --recursive -e robots=off --reject "index.html" --no-host-directories --cut-dirs=6 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6

grep -v '^#' GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_assembly_report.txt | cut -f7,10 | sed 's/\r$//' > chromosome.UCSC.txt

# Generate Annomen table
perl Annomen_table.pl -c chromosome.UCSC.txt GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz genome.fasta refseq.transcript.fasta refseq.protein.fasta refseq.transcript.gbff > Annomen_table.txt 2> Annomen_table.log
