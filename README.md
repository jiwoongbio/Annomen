# Annomen

Annotate variant nomenclature


## Annotations

* Gene name
* Region type
  * CDS
  * utr5
  * utr3
  * intron
  * non_coding_exon
  * non_coding_intron
* Mutation class
  * silent
  * missense
  * inframe
  * frameshift
  * nonsense
  * readthrough
  * startcodon
  * splicing
  * junction
* Strand
* Splice distance
* Transcript ID
* Protein ID
* Transcript variation nomenclature
* Protein variation nomenclature


## Requirements

1. Perl: https://www.perl.org
2. BioPerl: http://www.bioperl.org/wiki/Main_Page
   * Bio::DB::Fasta
   * Bio::SeqIO
3. EMBOSS: http://emboss.sourceforge.net
   * needle
   * stretcher
4. Basic linux commands: bash, gzip, sort, awk, sed, ...
5. lftp: http://lftp.yar.ru


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.
```
git clone https://github.com/jiwoongbio/Annomen.git
```


## Usages

1. Prepare annotation table file
   1. Download GTF file from http://genome.ucsc.edu/cgi-bin/hgTables
      * clade: Mammal
      * genome: Human
      * assembly: Feb. 2009 (GRCh37/hg19)
      * group: Genes and Gene Predictions
      * track: NCBI RefSeq
      * table: RefSeq All (ncbiRefSeq)
      * region: genome
      * output format: GTF - gene transfer format (limited)
      * output file: hg19_RefSeq.gtf.gz
      * file type returned: gzip compressed
      * get output
   2. Execute **Annomen_table.hg19.sh**
      ```
      ./Annomen_table.hg19.sh
      ```

2. Annotate variant file in VCF or tab-separated columns of chromosome, position, reference base, variant base
   * Execute **Annomen.hg19.sh**
     ```
     ./Annomen.hg19.sh <input file>
     ```
   * **Example**: annotating ClinVar variants
     ```
     ./clinvar.sh
     ```
