#!/usr/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

# Requirements
# 1. Perl: https://www.perl.org
# 2. BioPerl: http://www.bioperl.org/wiki/Main_Page
#    - Bio::DB::Fasta
#    - Bio::SeqIO
# 3. EMBOSS: http://emboss.sourceforge.net
#    - needle
#    - stretcher
# 4. Basic linux commands: bash, rm, gzip, sort, echo, find, sed, awk
# 5. lftp: http://lftp.yar.ru

# Remove old files
rm -rf hg38.fa.gz refseq/H_sapiens/mRNA_Prot human.rna.fna human.protein.faa human.rna.gbff Annomen_table.hg38.txt Annomen_table.hg38.log

# Prepare reference genome fasta file
genomeFastaFile=`readlink -f ~/data/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa`
if test -r $genomeFastaFile; then
	ln -sf $genomeFastaFile hg38.fasta
else 
	lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
	gzip -dc hg38.fa.gz > hg38.fasta
fi

# Download RefSeq files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot refseq/H_sapiens/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz;     do gzip -dc $file; done > human.rna.fna
for file in refseq/H_sapiens/mRNA_Prot/human.*.protein.faa.gz; do gzip -dc $file; done > human.protein.faa
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - human.rna.gbff

lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13'

grep -v '^#' GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt | cut -f7,10 | sed 's/\r$//' > chromosome.hg38.txt

# Generate Annomen table
perl Annomen_table.pl -c chromosome.hg38.txt GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz hg38.fasta human.rna.fna human.protein.faa human.rna.gbff > Annomen_table.hg38.txt 2> Annomen_table.hg38.log
