# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

# Requirements
# 1. Perl: https://www.perl.org
# 2. BioPerl: http://www.bioperl.org/wiki/Main_Page
#    - Bio::DB::Fasta
#    - Bio::SeqIO
# 3. EMBOSS: http://emboss.sourceforge.net
#    - needle
#    - stretcher
# 4. Basic linux commands: rm, gzip, awk, sort, echo, find
# 5. lftp: http://lftp.yar.ru

# http://genome.ucsc.edu/cgi-bin/hgTables
# - clade: Mammal
# - genome: Rat
# - assembly: Jul. 2014 (RGSC 6.0/rn6)
# - group: Genes and Gene Predictions
# - track: NCBI RefSeq
# - table: RefSeq All (ncbiRefSeq)
# - region: genome
# - output format: GTF - gene transfer format (limited)
# - output file: rn6_RefSeq.gtf.gz
# - file type returned: gzip compressed
# - get output

if [ ! -f "rn6_RefSeq.gtf.gz" ]; then
	echo "rn6_RefSeq.gtf.gz is not available." >&2
	exit 1
fi

# Remove old files
rm -rf rn6.fa.gz gene_info.gz gene2refseq.gz rn6_RefSeq.gtf refseq/R_norvegicus/mRNA_Prot rat.rna.fna rat.protein.faa rat.rna.gbff Annomen_table.txt Annomen_table.log

# Prepare reference genome fasta file
lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz'
gzip -dc rn6.fa.gz > rn6.fasta

# Download gene files from NCBI FTP
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz'

# Generate GTF file with gene names
perl RefSeq.gtf.pl rn6_RefSeq.gtf.gz gene_info.gz gene2refseq.gz 10116 > rn6_RefSeq.gtf

# Download RefSeq files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/R_norvegicus/mRNA_Prot refseq/R_norvegicus/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/R_norvegicus/mRNA_Prot/rat.*.rna.fna.gz;     do gzip -dc $file; done > rat.rna.fna
for file in refseq/R_norvegicus/mRNA_Prot/rat.*.protein.faa.gz; do gzip -dc $file; done > rat.protein.faa
for file in refseq/R_norvegicus/mRNA_Prot/rat.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - rat.rna.gbff

# Generate Annomen table
perl Annomen_table.pl rn6_RefSeq.gtf rn6.fasta rat.rna.fna rat.protein.faa rat.rna.gbff > Annomen_table.rn6.txt 2> Annomen_table.rn6.log
