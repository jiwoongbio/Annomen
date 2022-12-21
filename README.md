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
3. EMBOSS: http://emboss.sourceforge.net or EMBOSS-6.6.0.reduced.tar.gz
   * needle
   * stretcher
4. Basic linux commands: bash, rm, gzip, sort, echo, find, sed, awk
5. lftp: http://lftp.yar.ru


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.
```
git clone https://github.com/jiwoongbio/Annomen.git
```


## Usages

1. Prepare annotation table file
   * Execute **Annomen_table.hg38.sh**
      ```
      ./Annomen_table.hg38.sh
      ```

2. Annotate variant file in VCF or tab-separated columns of chromosome, position, reference base, variant base
   * Execute **Annomen.hg38.sh**
     ```
     ./Annomen.hg38.sh <input file>
     ```
   * **Example**: annotating ClinVar variants
     ```
     ./clinvar.sh
     ```
