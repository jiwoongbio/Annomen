# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

directory=`dirname $0`
inputFile=$1
referenceFastaFile=$directory/hg38.fasta
annotationTableFile=$directory/Annomen_table.hg38.txt
transcriptFastaFile=$directory/human.rna.fna
proteinFastaFile=$directory/human.protein.faa

perl $directory/leftalignIndel.pl $inputFile $referenceFastaFile | perl $directory/sort_by_reference.pl -c - $referenceFastaFile 0 1 | perl $directory/Annomen.pl - $referenceFastaFile $annotationTableFile $transcriptFastaFile $proteinFastaFile
