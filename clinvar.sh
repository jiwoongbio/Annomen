#!/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

directory=`dirname $0`

# Download VCF file from NCBI FTP and convert chromosome names
wget --no-verbose --no-check-certificate https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gzip -dc clinvar.vcf.gz | awk -F'\t' '($o ~ /^#/ || $1 ~ /^[0-9]+|X|Y|MT$/)' | sed 's/^/chr/' | sed 's/^chr#/#/' | sed 's/^chrMT/chrM/' > clinvar.vcf

# Annotate variant nomenclature
$directory/Annomen.sh clinvar.vcf > clinvar.annotated.vcf

# Generate variant table
perl $directory/vcf.table.pl -c column.name.txt -P clinvar.annotated.vcf > clinvar.annotated.txt
