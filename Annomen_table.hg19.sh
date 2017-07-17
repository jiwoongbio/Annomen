# Author: Jiwoong Kim (jiwoongbio@gmail.com)

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
# - assembly: Feb. 2009 (GRCh37/hg19)
# - group: Genes and Gene Predictions
# - track: RefSeq Genes
# - table: refGene
# - region: genome
# - output format: GTF - gene transfer format
# - output file: hg19_refGene.gtf.gz
# - file type returned: gzip compressed
# - get output

if [ ! -f "hg19_refGene.gtf.gz" ]; then
	echo "hg19_refGene.gtf.gz is not available."
	exit 1
fi

# Remove old files
rm -rf hg19_chromFa.tar.gz hg19_chromFa gene_info.gz gene2refseq.gz hg19_refGene.gtf refseq/H_sapiens/mRNA_Prot human.rna.fna human.protein.faa human.rna.gbff Annomen_table.txt Annomen_table.log

# Prepare reference genome fasta file
lftp -c 'get http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -o hg19_chromFa.tar.gz'
mkdir hg19_chromFa
tar -zxf hg19_chromFa.tar.gz -C hg19_chromFa
cat hg19_chromFa/chrM.fa `ls hg19_chromFa/chr*.fa | awk '($o ~ /^hg19_chromFa\/chr[0-9]+\.fa$/)' | sed 's/^hg19_chromFa\/chr//' | sort -n | sed 's/^/hg19_chromFa\/chr/'` hg19_chromFa/chrX.fa hg19_chromFa/chrY.fa hg19_chromFa/chrM_*.fa `ls hg19_chromFa/chr*.fa | awk '($o ~ /^hg19_chromFa\/chr[0-9]+_/)' | sed 's/^hg19_chromFa\/chr//' | sort -n | sed 's/^/hg19_chromFa\/chr/'` hg19_chromFa/chrX_*.fa hg19_chromFa/chrY_*.fa `ls hg19_chromFa/chr*.fa | awk '($o !~ /^hg19_chromFa\/chr(M|[0-9]+|X|Y)\.fa$/ && $o !~ /^hg19_chromFa\/chr(M|[0-9]+|X|Y)_/)'` > hg19.fasta

# Download gene files from NCBI FTP
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
lftp -c 'get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz'

# Generate GTF file with gene names
perl refGene.gtf.pl hg19_refGene.gtf.gz gene_info.gz gene2refseq.gz 9606 > hg19_refGene.gtf

# Download refGene files from NCBI FTP
lftp -c 'mirror -p -L ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot refseq/H_sapiens/mRNA_Prot'

# Prepare input files: 1. transcript FASTA, 2. protein FASTA, 3. transcript GenBank
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz;     do gzip -dc $file; done | sed 's/^>ref|/>/' | sed 's/|.*$//' | tr '\n' ' ' | sed 's/ $/\n/' | sed 's/ >/\n>/g' | grep -v '^>XM_' | grep -v '^>XR_' | sed 's/ /\n/g' > human.rna.fna
for file in refseq/H_sapiens/mRNA_Prot/human.*.protein.faa.gz; do gzip -dc $file; done | sed 's/^>ref|/>/' | sed 's/|.*$//' | tr '\n' ' ' | sed 's/ $/\n/' | sed 's/ >/\n>/g' | grep -v '^>XP_' | grep -v '^>YP_' | sed 's/ /\n/g' > human.protein.faa
for file in refseq/H_sapiens/mRNA_Prot/human.*.rna.gbff.gz;    do gzip -dc $file; done | perl splitGenBank.pl - human.rna.gbff

# Generate Annomen table
perl Annomen_table.pl hg19_refGene.gtf hg19.fasta human.rna.fna human.protein.faa human.rna.gbff > Annomen_table.txt 2> Annomen_table.log
