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
# - genome: Mouse
# - assembly: Dec. 2011 (GRCm38/mm10)
# - group: Genes and Gene Predictions
# - track: NCBI RefSeq
# - table: RefSeq All (ncbiRefSeq)
# - region: genome
# - output format: GTF - gene transfer format (limited)
# - output file: mm10_RefSeq.gtf.gz
# - file type returned: gzip compressed
# - get output

if [ ! -f "mm10_RefSeq.gtf.gz" ]; then
	echo "mm10_RefSeq.gtf.gz is not available." >&2
	exit 1
fi

# Remove old files
rm -rf mm10_chromFa.tar.gz mm10_chromFa gene_info.gz gene2refseq.gz mm10_RefSeq.gtf refseq/M_musculus/mRNA_Prot mouse.rna.fna mouse.protein.faa mouse.rna.gbff Annomen_table.txt Annomen_table.log

# Prepare reference genome fasta file
lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz -o mm10_chromFa.tar.gz'
mkdir mm10_chromFa
tar -zxf mm10_chromFa.tar.gz -C mm10_chromFa
cat mm10_chromFa/chrM.fa `ls mm10_chromFa/chr*.fa | awk '($o ~ /^mm10_chromFa\/chr[0-9]+\.fa$/)' | sed 's/^mm10_chromFa\/chr//' | sort -n | sed 's/^/mm10_chromFa\/chr/'` mm10_chromFa/chrX.fa mm10_chromFa/chrY.fa mm10_chromFa/chrM_*.fa `ls mm10_chromFa/chr*.fa | awk '($o ~ /^mm10_chromFa\/chr[0-9]+_/)' | sed 's/^mm10_chromFa\/chr//' | sort -n | sed 's/^/mm10_chromFa\/chr/'` mm10_chromFa/chrX_*.fa mm10_chromFa/chrY_*.fa `ls mm10_chromFa/chr*.fa | awk '($o !~ /^mm10_chromFa\/chr(M|[0-9]+|X|Y)\.fa$/ && $o !~ /^mm10_chromFa\/chr(M|[0-9]+|X|Y)_/)'` > mm10.fasta

# Download gene files from NCBI FTP
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz'

# Generate GTF file with gene names
perl RefSeq.gtf.pl mm10_RefSeq.gtf.gz gene_info.gz gene2refseq.gz 10090 > mm10_RefSeq.gtf

# Download RefSeq files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot refseq/M_musculus/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.fna.gz;     do gzip -dc $file; done > mouse.rna.fna
for file in refseq/M_musculus/mRNA_Prot/mouse.*.protein.faa.gz; do gzip -dc $file; done > mouse.protein.faa
for file in refseq/M_musculus/mRNA_Prot/mouse.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - mouse.rna.gbff

# Generate Annomen table
perl Annomen_table.pl mm10_RefSeq.gtf mm10.fasta mouse.rna.fna mouse.protein.faa mouse.rna.gbff > Annomen_table.mm10.txt 2> Annomen_table.mm10.log
