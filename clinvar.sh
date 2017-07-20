# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

directory=`dirname $0`

# Download VCF file from NCBI FTP and convert chromosome names
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
gzip -dc clinvar.vcf.gz | awk -F'\t' '($o ~ /^#/ || $1 ~ /^[0-9]+|X|Y$/)' | sed 's/^/chr/' | sed 's/^chr#/#/' > clinvar.vcf

# Annotate variant nomenclature
$directory/Annomen.hg19.sh clinvar.vcf > clinvar.annotated.vcf

# Generate variant table
perl $directory/vcf.table.pl -p clinvar.annotated.vcf `grep -v '^#' column.txt | cut -f2` | bash -c "cat <(grep -v '^#' column.txt | cut -f1 | tr '\n' '\t' | sed 's/\t$/\n/' | sed 's/^/#/') -" > clinvar.annotated.txt
