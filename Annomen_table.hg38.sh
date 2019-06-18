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
# - genome: Human
# - assembly: Dec. 2013 (GRCh38/hg38)
# - group: Genes and Gene Predictions
# - track: NCBI RefSeq
# - table: RefSeq All (ncbiRefSeq)
# - region: genome
# - output format: GTF - gene transfer format (limited)
# - output file: hg38_RefSeq.gtf.gz
# - file type returned: gzip compressed
# - get output

if [ ! -f "hg38_RefSeq.gtf.gz" ]; then
	echo "hg38_RefSeq.gtf.gz is not available." >&2
	exit 1
fi

# Remove old files
rm -rf hg38.fa.gz gene_info.gz gene2refseq.gz hg38_RefSeq.gtf refseq/H_sapiens/mRNA_Prot human.rna.fna human.protein.faa human.rna.gbff Annomen_table.txt Annomen_table.log

# Prepare reference genome fasta file
lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
gzip -dc hg38.fa.gz > hg38.fasta

# Download gene files from NCBI FTP
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz'

# Generate GTF file with gene names
perl RefSeq.gtf.pl hg38_RefSeq.gtf.gz gene_info.gz gene2refseq.gz 9606 > hg38_RefSeq.gtf

# Download RefSeq files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot refseq/H_sapiens/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz;     do gzip -dc $file; done > human.rna.fna
for file in refseq/H_sapiens/mRNA_Prot/human.*.protein.faa.gz; do gzip -dc $file; done > human.protein.faa
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - human.rna.gbff

# Generate Annomen table
perl Annomen_table.pl hg38_RefSeq.gtf hg38.fasta human.rna.fna human.protein.faa human.rna.gbff > Annomen_table.hg38.txt 2> Annomen_table.hg38.log
