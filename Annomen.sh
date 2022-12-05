#!/usr/bin/bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

directory=`dirname $0`
inputFile=$1
referenceFastaFile=$directory/genome.fasta
annotationTableFile=$directory/Annomen_table.txt
transcriptFastaFile=$directory/refseq.transcript.fasta
proteinFastaFile=$directory/refseq.protein.fasta

perl $directory/leftalignIndel.pl $inputFile $referenceFastaFile | perl $directory/sort_by_reference.pl - $referenceFastaFile 0 1 | perl $directory/Annomen.pl - $referenceFastaFile $annotationTableFile $transcriptFastaFile $proteinFastaFile
