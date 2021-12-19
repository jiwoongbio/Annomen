#!/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

directory=`dirname $0`

# Download VCF file from NCBI FTP and convert chromosome names
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
gzip -dc clinvar.vcf.gz | awk -F'\t' '($o ~ /^#/ || $1 ~ /^[0-9]+|X|Y$/)' | sed 's/^/chr/' | sed 's/^chr#/#/' > clinvar.vcf

# Annotate variant nomenclature
$directory/Annomen.hg38.sh clinvar.vcf > clinvar.annotated.vcf

# Generate variant table
perl $directory/vcf.table.pl -c column.name.txt -P clinvar.annotated.vcf > clinvar.annotated.txt
