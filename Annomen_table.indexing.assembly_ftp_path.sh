#!/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

# Requirements
# 1. Perl: https://www.perl.org
# 2. BioPerl: http://www.bioperl.org/wiki/Main_Page
#    - Bio::DB::Fasta
#    - Bio::SeqIO
# 3. EMBOSS: http://emboss.sourceforge.net or EMBOSS-6.6.0.reduced.tar.gz
#    - needle
#    - stretcher
# 4. Basic linux commands: bash, rm, gzip, sort, echo, find, sed, awk, wget

assembly_ftp_path=$1
assembly_prefix=$(echo $assembly_ftp_path | sed 's/^.*\///')

wget --no-verbose --no-check-certificate --recursive --execute robots=off --reject 'index.html*' --no-parent --level=1 --no-host-directories --cut-dirs=6 $assembly_ftp_path/

(file=${assembly_prefix}/${assembly_prefix}_genomic.fna.gz; [ -e $file ] && gzip -dc $file > genome.fasta)
(file=${assembly_prefix}/${assembly_prefix}_rna.fna.gz; [ -e $file ] && gzip -dc $file > transcript.fasta)
(file=${assembly_prefix}/${assembly_prefix}_protein.faa.gz; [ -e $file ] && gzip -dc $file > protein.fasta)

# Remove old Annomen table files
rm -rf Annomen_table.txt Annomen_table.log

# Generate Annomen table files
perl Annomen_table.pl ${assembly_prefix}/${assembly_prefix}_genomic.gff.gz genome.fasta transcript.fasta protein.fasta > Annomen_table.txt 2> Annomen_table.log


# Indexing

# Generate genome index
time bwa index genome.fasta
time samtools faidx genome.fasta

#PICARD=where/is/picard.jar
time java -jar $PICARD CreateSequenceDictionary REFERENCE=genome.fasta

# Calculate genome sequence lengths
time fasta.length.pl genome.fasta > genome.length.txt

# Generate transcript index
time bwa index transcript.fasta
time samtools faidx transcript.fasta

# Calculate transcript sequence lengths
time fasta.length.pl transcript.fasta > transcript.length.txt

# Unzip GTF file
time gzip -dc */*_genomic.gtf.gz | grep -v '^#' | sed -r 's/""([^;][^"]*)"";/"\1";/' > genome.gtf

# Generate STAR index
rm -rf STAR; mkdir STAR; time STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles genome.fasta --sjdbGTFfile genome.gtf

# List rRNA
time perl gff_extract.pl -E -f feature=rRNA -f gene_biotype=rRNA -b or */*_genomic.gff.gz transcript_id,locus_tag | sort -u > rRNA.txt

# List locus tag / product description pairs
time perl gff_extract.pl -E */*_genomic.gff.gz transcript_id,locus_tag product | sort -u > product.txt
