# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

wget ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-91/variation/vcf/mus_musculus/mus_musculus.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz

gzip -d Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gzip -d mus_musculus.vcf.gz
gzip -d Mus_musculus.GRCm38.91.gtf.gz
gzip -d Mus_musculus.GRCm38.cdna.all.fa.gz
gzip -d Mus_musculus.GRCm38.pep.all.fa.gz

sed 's/\.[0-9]* .*$//' Mus_musculus.GRCm38.cdna.all.fa > Mus_musculus.GRCm38.cdna.all.header_modified.fa
sed 's/\.[0-9]* .*$//' Mus_musculus.GRCm38.pep.all.fa > Mus_musculus.GRCm38.pep.all.header_modified.fa
perl Annomen_table.pl Mus_musculus.GRCm38.91.gtf Mus_musculus.GRCm38.dna.primary_assembly.fa Mus_musculus.GRCm38.cdna.all.header_modified.fa Mus_musculus.GRCm38.pep.all.header_modified.fa > Annomen_table.txt 2> Annomen_table.log

perl leftalignIndel.pl mus_musculus.vcf Mus_musculus.GRCm38.dna.primary_assembly.fa | perl sort_by_reference.pl -c - Mus_musculus.GRCm38.dna.primary_assembly.fa 0 1 | grep -v '^#' | awk -F'\t' -vOFS='\t' 'BEGIN {print "chromosome", "start", "end", "haplotypeReference", "haplotypeAlternate", "name"} {print $1, $2, $2 + length($4) - 1, $4, $5, $3}' > mus_musculus.variation.txt
