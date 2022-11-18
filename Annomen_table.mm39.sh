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
rm -rf mm39.fa.gz refseq/M_musculus/mRNA_Prot mouse.rna.fna mouse.protein.faa mouse.rna.gbff Annomen_table.mm39.txt Annomen_table.mm39.log

# Prepare reference genome fasta file
if test -r "$genomeFastaFile"; then
	ln -sf "$genomeFastaFile" mm39.fasta
else 
	lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz'
	gzip -d mm39.fa.gz
	ln -sf mm39.fa mm39.fasta
fi

# Download RefSeq files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot refseq/M_musculus/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.fna.gz;     do gzip -dc $file; done > mouse.rna.fna
for file in refseq/M_musculus/mRNA_Prot/mouse.*.protein.faa.gz; do gzip -dc $file; done > mouse.protein.faa
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - mouse.rna.gbff

lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39'

grep -v '^#' GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt | cut -f7,10 | sed 's/\r$//' > chromosome.mm39.txt

# Generate Annomen table
perl Annomen_table.pl -c chromosome.mm39.txt GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz mm39.fasta mouse.rna.fna mouse.protein.faa mouse.rna.gbff > Annomen_table.mm39.txt 2> Annomen_table.mm39.log
