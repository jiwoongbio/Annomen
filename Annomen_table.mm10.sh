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
rm -rf mm10.fa.gz refseq/M_musculus/mRNA_Prot mouse.rna.fna mouse.protein.faa mouse.rna.gbff Annomen_table.mm10.txt Annomen_table.mm10.log

# Prepare reference genome fasta file
if test -r "$genomeFastaFile"; then
	ln -sf $genomeFastaFile mm10.fasta
else 
	lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz'
	gzip -d mm10.fa.gz
	ln -sf mm10.fa mm10.fasta
fi

# Download RefSeq files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot refseq/M_musculus/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.fna.gz;     do gzip -dc $file; done > mouse.rna.fna
for file in refseq/M_musculus/mRNA_Prot/mouse.*.protein.faa.gz; do gzip -dc $file; done > mouse.protein.faa
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - mouse.rna.gbff

lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6'

grep -v '^#' GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_assembly_report.txt | cut -f7,10 | sed 's/\r$//' > chromosome.mm10.txt

# Generate Annomen table
perl Annomen_table.pl -c chromosome.mm10.txt GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz mm10.fasta mouse.rna.fna mouse.protein.faa mouse.rna.gbff > Annomen_table.mm10.txt 2> Annomen_table.mm10.log
