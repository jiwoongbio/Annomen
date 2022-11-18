#!/usr/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

directory=`dirname $0`
inputFile=$1
referenceFastaFile=$directory/mm39.fasta
annotationTableFile=$directory/Annomen_table.mm39.txt
transcriptFastaFile=$directory/mouse.rna.fna
proteinFastaFile=$directory/mouse.protein.faa

perl $directory/leftalignIndel.pl $inputFile $referenceFastaFile | perl $directory/sort_by_reference.pl - $referenceFastaFile 0 1 | perl $directory/Annomen.pl - $referenceFastaFile $annotationTableFile $transcriptFastaFile $proteinFastaFile
