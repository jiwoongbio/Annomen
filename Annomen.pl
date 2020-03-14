# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum min max);
use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	'C=s' => \@codonList,
	'S=s' => \(my $startCodons = 'ATG,CTG,TTG'),
	'r' => \(my $checkReferenceVariation = ''),
	'p=i' => \(my $peptideFlankingLength = 0),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl Annomen.pl [options] variant.vcf|variant.txt genome.fasta Annomen_table.txt transcript.fasta protein.fasta > variant.annotated.vcf|variant.annotated.txt

Options: -h       display this help message
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 1 (standard)]
         -S STR   comma-separated start codons [$startCodons]
         -r       check reference variation
         -p       peptide flanking length

EOF
}
{
	my %codonHash = (
		'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
		'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
		'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
		'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

		'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
		'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
		'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
		'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

		'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
		'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
		'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
		'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

		'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
		'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
		'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
		'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
	);
	$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);

	sub translate {
		my ($sequence) = @_;
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. int(length($sequence) / 3) - 1);
	}
}
my %startCodonHash = map {$_ => 1} split(/,/, $startCodons);
my ($variantFile, $referenceFastaFile, $annotationTableFile, $transcriptFastaFile, $proteinFastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
my %aaOneToThreeLetter = ('A' => 'Ala', 'B' => 'Asx', 'C' => 'Cys', 'D' => 'Asp', 'E' => 'Glu', 'F' => 'Phe', 'G' => 'Gly', 'H' => 'His', 'I' => 'Ile', 'K' => 'Lys', 'L' => 'Leu', 'M' => 'Met', 'N' => 'Asn', 'P' => 'Pro', 'Q' => 'Gln', 'R' => 'Arg', 'S' => 'Ser', 'T' => 'Thr', 'V' => 'Val', 'W' => 'Trp', 'X' => 'Xxx', 'Y' => 'Tyr', 'Z' => 'Glx');
my %transcriptSequenceHash = ();
{
	my $transcriptId = '';
	open(my $reader, ($transcriptFastaFile =~ /\.gz$/) ? "gzip -dc $transcriptFastaFile |" : $transcriptFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($transcriptId = $1));
		$transcriptSequenceHash{$transcriptId} .= $line;
	}
	close($reader);
}
my %proteinSequenceHash = ();
{
	my $proteinId = '';
	open(my $reader, ($proteinFastaFile =~ /\.gz$/) ? "gzip -dc $proteinFastaFile |" : $proteinFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($proteinId = $1));
		$proteinSequenceHash{$proteinId} .= $line;
	}
	close($reader);
	s/U/*/g foreach(values %proteinSequenceHash);
	s/\*?$/*/ foreach(values %proteinSequenceHash);
}
{
	my @chromosomeList = getChromosomeList();
	my %chromosomeIndexHash = ();
	@chromosomeIndexHash{@chromosomeList} = 0 .. scalar(@chromosomeList) - 1;
	open(my $reader, $annotationTableFile);
	sub getTokenList {
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^#/);
			return split(/\t/, $line, -1);
		}
		return ();
	}
	my @columnList = getTokenList;
	my %columnIndexHash = ();
	@columnIndexHash{@columnList} = 0 .. $#columnList;
	my @tokenList = getTokenList;
	my @tokenListList = ();
	sub getTokenListList {
		my ($chromosome, $start, $end) = @_;
		@tokenListList = grep {$_->[$columnIndexHash{'chromosome'}] eq $chromosome && $start <= $_->[$columnIndexHash{'end'}]} @tokenListList;
		while(@tokenList && defined($chromosomeIndexHash{$chromosome}) && $chromosomeIndexHash{$tokenList[$columnIndexHash{'chromosome'}]} < $chromosomeIndexHash{$chromosome}) {
			@tokenList = getTokenList;
		}
		while(@tokenList && $tokenList[$columnIndexHash{'chromosome'}] eq $chromosome && $tokenList[$columnIndexHash{'start'}] <= $end) {
			push(@tokenListList, [@tokenList]) if($start <= $tokenList[$columnIndexHash{'end'}]);
			@tokenList = getTokenList;
		}
		return grep {$_->[$columnIndexHash{'start'}] <= $end} @tokenListList;
	}
	sub getColumnList {
		return @columnList;
	}
	sub closeAnnotationTableFileReader {
		close($reader);
	}
}
my $title = 'Annomen';
my @columnList = ('allele', 'geneName', 'region', 'mutation', 'strand', 'spliceDistance', 'transcriptId', 'proteinId', 'nucleotideVariationNomenclature', 'proteinVariationNomenclature', ($peptideFlankingLength ? ('peptidePosition', 'originalPeptide', 'mutationPeptide') : ()));
my $vcf = '';
open(my $reader, $variantFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^#/) {
		$vcf = 1 if($line =~ /^##fileformat=VCF/);
		print "$line\n";
		next;
	}
	my @tokenList = split(/\t/, $line, -1);
	my ($chromosome, $position, $refBase, $altBase) = @tokenList[$vcf ? (0, 1, 3, 4) : (0, 1, 2, 3)];
	($refBase, $altBase) = map {uc} ($refBase, $altBase);
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	my @tokenHashList = ();
	foreach my $altBaseIndex (0 .. $#altBaseList) {
		if((my $altBase = $altBaseList[$altBaseIndex]) =~ /^[A-Za-z]*$/) {
			my ($start, $end) = getExtendedStartEnd($chromosome, stripIdentical($position, $refBase, $altBase));
			my @transcriptIdList = ();
			my %transcriptIdTokenHashListHash = ();
			foreach my $tokenList (getTokenListList($chromosome, $position, $end)) {
				my %tokenHash = ();
				$tokenHash{'allele'} = $altBaseIndex + 1 if(scalar(@altBaseList) > 1);
				@tokenHash{(getColumnList)} = @$tokenList;
				next if($tokenHash{'end'} < $start);
				if(($tokenHash{'spliceDistance'} = getSpliceDistance($start, $end, \%tokenHash)) != 0) {
					setMutationNomenclatures($chromosome, $position, $refBase, $altBase, \%tokenHash);
				}
				$tokenHash{'region'} = getRegion($tokenHash{'strand'}, $start, $end, @tokenHash{'region', 'codingStart', 'codingEnd'});
				my $transcriptId = $tokenHash{'transcriptId'};
				push(@transcriptIdList, $transcriptId) unless(defined($transcriptIdTokenHashListHash{$transcriptId}));
				push(@{$transcriptIdTokenHashListHash{$transcriptId}}, \%tokenHash);
			}
			foreach my $transcriptId (@transcriptIdList) {
				my %tokenHash = %{$transcriptIdTokenHashListHash{$transcriptId}->[0]};
				if(scalar(my @tokenHashList = @{$transcriptIdTokenHashListHash{$transcriptId}}) > 1) {
					my @regionList = map {split(/\^/, $_->{'region'})} sort {$a->{'startInTranscript'} <=> $b->{'startInTranscript'}} @tokenHashList;
					$tokenHash{'region'} = join('^', map {$regionList[$_]} grep {$_ == 0 || $regionList[$_ - 1] ne $regionList[$_]} 0 .. $#regionList);
					$tokenHash{'spliceDistance'} = min(map {$_->{'spliceDistance'}} @tokenHashList);
				}
				$tokenHash{'mutation'} = 'junction' if($tokenHash{'spliceDistance'} == 0);
				push(@tokenHashList, \%tokenHash);
			}
		}
	}
	if($vcf) {
		if(@tokenHashList) {
			my @infoList = grep {$_ ne '.' && $_ !~ /^$title\./} split(/;/, $tokenList[7]);
			for(my $number = 1; $number <= scalar(@tokenHashList); $number++) {
				my %tokenHash = %{$tokenHashList[$number - 1]};
				if(scalar(@tokenHashList) > 1) {
					push(@infoList, map {"$title.${_}_$number=$tokenHash{$_}"} grep {defined($tokenHash{$_}) && $tokenHash{$_} ne ''} @columnList);
				} else {
					push(@infoList, map {"$title.$_=$tokenHash{$_}"} grep {defined($tokenHash{$_}) && $tokenHash{$_} ne ''} @columnList);
				}
			}
			$tokenList[7] = join(';', @infoList);
		}
		print join("\t", @tokenList), "\n";
	} else {
		push(@tokenHashList, {}) if(scalar(@tokenHashList) == 0);
		open(my $writer, "| sort -u");
		print $writer join("\t", @tokenList, map {defined($_) ? $_ : ''} @$_{@columnList}), "\n" foreach(@tokenHashList);
		close($writer);
	}
}
close($reader);
closeAnnotationTableFileReader;

sub getRegion {
	my ($strand, $start, $end, $region, $codingStart, $codingEnd) = @_;
	if(grep {$_ ne ''} ($codingStart, $codingEnd)) {
		if($region eq 'exon') {
			if($codingStart <= $start && $end <= $codingEnd) {
				return 'CDS';
			} elsif($end < $codingStart) {
				return 'utr5' if($strand eq '+');
				return 'utr3' if($strand eq '-');
			} elsif($codingEnd < $start) {
				return 'utr3' if($strand eq '+');
				return 'utr5' if($strand eq '-');
			} elsif($start < $codingStart && $codingEnd < $end) {
				return 'utr5^CDS^utr3';
			} elsif($start < $codingStart) {
				return 'utr5^CDS' if($strand eq '+');
				return 'CDS^utr3' if($strand eq '-');
			} elsif($codingEnd < $end) {
				return 'CDS^utr3' if($strand eq '+');
				return 'utr5^CDS' if($strand eq '-');
			}
		}
	} else {
		return "non_coding_$region";
	}
	return $region;
}

sub getSpliceDistance {
	my ($start, $end, $tokenHash) = @_;
	my $forwardSpliceDistance = $start - $tokenHash->{'start'} + 1;
	my $reverseSpliceDistance = $tokenHash->{'end'} - $end + 1;
	$forwardSpliceDistance = 0 if($forwardSpliceDistance < 0);
	$reverseSpliceDistance = 0 if($reverseSpliceDistance < 0);
	if($forwardSpliceDistance < $reverseSpliceDistance) {
		if($tokenHash->{'strand'} eq '+') {
			return $forwardSpliceDistance;
		} else {
			return -$forwardSpliceDistance;
		}
	}
	if($reverseSpliceDistance < $forwardSpliceDistance) {
		if($tokenHash->{'strand'} eq '-') {
			return $reverseSpliceDistance;
		} else {
			return -$reverseSpliceDistance;
		}
	}
	return $forwardSpliceDistance;
}

sub setMutationNomenclatures {
	my ($chromosome, $position, $refBase, $altBase, $tokenHash) = @_;
	($position, $refBase, $altBase) = stripIdentical($position, $refBase, $altBase);
	my @codingRegionList = map {[split(/\.\./, $_)]} split(/,/, $tokenHash->{'codingRegions'});
	@$tokenHash{'codingStartInTranscript', 'codingEndInTranscript'} = map {$codingRegionList[$_]->[$_]} (0, -1) if(@codingRegionList);
	if($tokenHash->{'region'} eq 'exon') {
		my @mismatchList = map {[split(':', $_, 3)]} split(/,/, $tokenHash->{'mismatch'});
		foreach my $mismatch (@mismatchList) {
			if($mismatch->[0] < $position && $position <= $mismatch->[0] + length($mismatch->[1])) {
				my $head = substr($mismatch->[1], 0, $position - $mismatch->[0]);
				($refBase, $altBase) = map {"$head$_"} ($refBase, $altBase);
				$position = $mismatch->[0];
			}
			if($mismatch->[0] <= $position + length($refBase) && $position + length($refBase) < $mismatch->[0] + length($mismatch->[1])) {
				my $tail = substr($mismatch->[1], $position + length($refBase) - ($mismatch->[0] + length($mismatch->[1])));
				($refBase, $altBase) = map {"$_$tail"} ($refBase, $altBase);
			}
		}
		my $rnaBase = $refBase;
		foreach my $mismatch (sort {$b->[0] <=> $a->[0]} @mismatchList) {
			if($position <= $mismatch->[0] && $mismatch->[0] + length($mismatch->[1]) <= $position + length($refBase)) {
				substr($rnaBase, $mismatch->[0] - $position, length($mismatch->[1]), $mismatch->[2]);
			}
		}
		$position = $position - $tokenHash->{'start'} + sum(0, map {length($_->[2]) - length($_->[1])} grep {$_->[0] < $position} @mismatchList);
		if($tokenHash->{'strand'} eq '+') {
			$position = $tokenHash->{'startInTranscript'} + $position;
		} else {
			$position = $tokenHash->{'endInTranscript'} - $position - length($rnaBase) + 1;
			($rnaBase, $refBase, $altBase) = map {getReverseComplementarySequence($_)} ($rnaBase, $refBase, $altBase);
		}
		my $transcriptSequence = $transcriptSequenceHash{my $transcriptId = $tokenHash->{'transcriptId'}};
		my @transcriptProteinVariationNomenclaturesList = ();
		push(@transcriptProteinVariationNomenclaturesList, [getTranscriptProteinVariationNomenclatures($transcriptSequence, $position, $rnaBase, $refBase, $tokenHash, @codingRegionList)]) if($checkReferenceVariation && $rnaBase ne $refBase);
		push(@transcriptProteinVariationNomenclaturesList, [getTranscriptProteinVariationNomenclatures($transcriptSequence, $position, $rnaBase, $altBase, $tokenHash, @codingRegionList)]);
		$tokenHash->{'nucleotideVariationNomenclature'} = join(',', map {$_->[0]} @transcriptProteinVariationNomenclaturesList);
		if(grep {$_ ne ''} (my @proteinVariationNomenclatureList = map {$_->[1]} @transcriptProteinVariationNomenclaturesList)) {
			$tokenHash->{'proteinVariationNomenclature'} = join(',', @proteinVariationNomenclatureList);
			$tokenHash->{'mutation'} = join(',', map {getProteinMutation($_)} @proteinVariationNomenclatureList);
		}
	}
	if($tokenHash->{'region'} eq 'intron') {
		if(abs($tokenHash->{'spliceDistance'}) <= 10) {
			my $intronSequence = uc($db->seq($chromosome, $tokenHash->{'start'}, $tokenHash->{'end'}));
			if($tokenHash->{'strand'} eq '+') {
				$position = $position - $tokenHash->{'start'} + 1;
			} else {
				$position = $tokenHash->{'end'} - $position + 1 - length($refBase) + 1;
				($intronSequence, $refBase, $altBase) = map {getReverseComplementarySequence($_)} ($intronSequence, $refBase, $altBase);
			}
			$tokenHash->{'intronSequenceLength'} = length($intronSequence);
#			if($tokenHash->{'spliceDistance'} > 0) {
#				$tokenHash->{'nucleotideVariationNomenclature'} = getNucleotideVariationNomenclature($intronSequence,  leftalignIndel($intronSequence, stripIdentical($position, $refBase, $altBase)), $tokenHash);
#			} else {
#				$tokenHash->{'nucleotideVariationNomenclature'} = getNucleotideVariationNomenclature($intronSequence, rightalignIndel($intronSequence, stripIdentical($position, $refBase, $altBase)), $tokenHash);
#			}
			$tokenHash->{'nucleotideVariationNomenclature'} = getNucleotideVariationNomenclature($intronSequence, rightalignIndel($intronSequence, stripIdentical($position, $refBase, $altBase)), $tokenHash);
		}
		$tokenHash->{'mutation'} = 'splicing' if(abs($tokenHash->{'spliceDistance'}) <= 2);
	}
}

sub getTranscriptProteinVariationNomenclatures {
	my ($transcriptSequence, $transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation, $tokenHash, @codingRegionList) = @_;
	my ($transcriptVariationNomenclature, $proteinVariationNomenclature) = ('', '');
	if($originalTranscriptVariation ne $mutationTranscriptVariation) {
#		($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation) = stripIdentical($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation);
#		if(@codingRegionList) {
#			($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation) =  leftalignIndel($transcriptSequence, $transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation);
#			($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation) = rightalignIndel($transcriptSequence, $transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation) if($tokenHash->{'codingStartInTranscript'} <= $transcriptVariationPosition);
#		} else {
#			($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation) = rightalignIndel($transcriptSequence, $transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation);
#		}
		($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation) = rightalignIndel($transcriptSequence, stripIdentical($transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation));
		$transcriptVariationNomenclature = getNucleotideVariationNomenclature($transcriptSequence, $transcriptVariationPosition, $originalTranscriptVariation, $mutationTranscriptVariation, $tokenHash);
		my ($transcriptVariationStartPosition, $transcriptVariationEndPosition) = sort {$a <=> $b} ($transcriptVariationPosition, $transcriptVariationPosition + length($originalTranscriptVariation) - 1);
		if(my ($index) = grep {$codingRegionList[$_]->[0] <= $transcriptVariationStartPosition && $transcriptVariationEndPosition <= $codingRegionList[$_]->[-1]} 0 .. $#codingRegionList) {
			my $codingVariationPosition = sum((map {$_->[-1] - $_->[0] + 1} @codingRegionList[0 .. $index - 1]), $transcriptVariationPosition - $codingRegionList[$index]->[0] + 1);
			my $originalCodingSequence = join('', (map {substr($transcriptSequence, $_->[0] - 1, $_->[-1] - $_->[0] + 1)} @codingRegionList), substr($transcriptSequence, $codingRegionList[-1]->[-1]));
			my $mutationCodingSequence = $originalCodingSequence;
			substr($mutationCodingSequence, $codingVariationPosition - 1, length($originalTranscriptVariation), $mutationTranscriptVariation);
			if(my $frame = $tokenHash->{'frame'}) {
				$codingVariationPosition = $codingVariationPosition - $frame;
				if($frame > 0) {
					$originalCodingSequence = substr($originalCodingSequence, $frame);
					$mutationCodingSequence = substr($mutationCodingSequence, $frame);
				}
				if($frame < 0) {
					$originalCodingSequence = ('N' x -$frame) . $originalCodingSequence;
					$mutationCodingSequence = ('N' x -$frame) . $mutationCodingSequence;
				}
			}
			my $proteinVariationPosition = int(($codingVariationPosition - 1) / 3) + 1;
			my $originalProteinSequence = $proteinSequenceHash{$tokenHash->{'proteinId'}};
			my $mutationProteinSequence;
			my ($originalProteinVariationLength, $mutationProteinVariationLength) = (0, 0);
			my $frameshift = (length($mutationCodingSequence) - length($originalCodingSequence)) % 3;
			if($frameshift == 0) {
				$originalProteinVariationLength = int(($codingVariationPosition + length($originalTranscriptVariation) - 2) / 3) + 2 - $proteinVariationPosition;
				$mutationProteinVariationLength = $originalProteinVariationLength + (length($mutationCodingSequence) - length($originalCodingSequence)) / 3;
				if($proteinVariationPosition - 1 + $originalProteinVariationLength <= length($originalProteinSequence)) {
					$mutationProteinSequence = $originalProteinSequence;
					my $originalProteinVariation = substr($originalProteinSequence, $proteinVariationPosition - 1, $originalProteinVariationLength);
					my $mutationProteinVariation = translate(substr($mutationCodingSequence, ($proteinVariationPosition - 1) * 3, $mutationProteinVariationLength * 3));
					if($mutationProteinVariation ne $originalProteinVariation) {
						($proteinVariationPosition, $originalProteinVariation, $mutationProteinVariation) = rightalignIndel($originalProteinSequence, stripIdentical($proteinVariationPosition, $originalProteinVariation, $mutationProteinVariation));
						$originalProteinVariationLength = length($originalProteinVariation);
						if(($mutationProteinVariationLength = index($mutationProteinVariation, '*') + 1) < 1) {
							$mutationProteinVariationLength = length($mutationProteinVariation);
						}
						substr($mutationProteinSequence, $proteinVariationPosition - 1, $originalProteinVariationLength, $mutationProteinVariation);
					}
				}
			}
			unless(defined($mutationProteinSequence)) {
				$mutationProteinSequence = substr($originalProteinSequence, 0, $proteinVariationPosition - 1) . translate(substr($mutationCodingSequence, ($proteinVariationPosition - 1) * 3));
			}
			if((my $index = index($originalProteinSequence, '*', $proteinVariationPosition - 1)) >= 0) {
				$originalProteinSequence = substr($originalProteinSequence, 0, $index);
			}
			if((my $index = index($mutationProteinSequence, '*', $proteinVariationPosition - 1)) >= 0) {
				$mutationProteinSequence = substr($mutationProteinSequence, 0, $index);
			}
			$originalProteinSequence = "$originalProteinSequence*";
			$mutationProteinSequence = "$mutationProteinSequence*";
			if($originalProteinSequence eq $mutationProteinSequence) {
				$proteinVariationNomenclature = 'p.=';
			} else {
				$proteinVariationPosition++ while(substr($originalProteinSequence, $proteinVariationPosition - 1, 1) eq substr($mutationProteinSequence, $proteinVariationPosition - 1, 1));
				my ($originalStartCodon, $mutationStartCodon) = map {substr($_, 0, 3)} ($originalCodingSequence, $mutationCodingSequence);
				if($originalStartCodon eq $mutationStartCodon || $startCodonHash{$mutationStartCodon}) {
					$proteinVariationNomenclature = getProteinVariationNomenclature($originalProteinSequence, $mutationProteinSequence, $proteinVariationPosition, $originalProteinVariationLength, $mutationProteinVariationLength, $frameshift);
					if($peptideFlankingLength) {
						my $index = ($proteinVariationPosition - 1) - $peptideFlankingLength;
						$index = 0 if($index < 0);
						my $originalPeptideLength = $frameshift == 0 ? $originalProteinVariationLength + $peptideFlankingLength * 2 : length($originalProteinSequence) - $index;
						my $mutationPeptideLength = $frameshift == 0 ? $mutationProteinVariationLength + $peptideFlankingLength * 2 : length($mutationProteinSequence) - $index;
						my $originalPeptide = substr($originalProteinSequence, $index, $originalPeptideLength);
						my $mutationPeptide = substr($mutationProteinSequence, $index, $mutationPeptideLength);
						$tokenHash->{'peptidePosition'} = $index + 1;
						$tokenHash->{'originalPeptide'} = $originalPeptide;
						$tokenHash->{'mutationPeptide'} = $mutationPeptide;
					}
				} else {
					$proteinVariationNomenclature = 'p.?';
				}
			}
		}
	}
	return ($transcriptVariationNomenclature, $proteinVariationNomenclature);
}

sub getProteinVariationNomenclature {
	my ($originalProteinSequence, $mutationProteinSequence, $proteinVariationPosition, $originalProteinVariationLength, $mutationProteinVariationLength, $frameshift) = @_;
	my @originalResidueList = map {$_ eq '*' ? $_ : $aaOneToThreeLetter{$_}} split(//, substr($originalProteinSequence, $proteinVariationPosition - 1));
	my @mutationResidueList = map {$_ eq '*' ? $_ : $aaOneToThreeLetter{$_}} split(//, substr($mutationProteinSequence, $proteinVariationPosition - 1));
	return "p.$originalResidueList[0]$proteinVariationPosition$mutationResidueList[0]" if($mutationResidueList[0] eq '*');
	return "p.$originalResidueList[0]$proteinVariationPosition$mutationResidueList[0]ext*" . (scalar(@mutationResidueList) - 1) if($originalResidueList[0] eq '*');
	return "p.$originalResidueList[0]$proteinVariationPosition$mutationResidueList[0]fs*" . scalar(@mutationResidueList) if($frameshift);
	my ($startNumber, $endNumber) = ($proteinVariationPosition, $proteinVariationPosition + $originalProteinVariationLength - 1);
	@originalResidueList = @originalResidueList[0 .. $originalProteinVariationLength - 1];
	@mutationResidueList = @mutationResidueList[0 .. $mutationProteinVariationLength - 1];
	return "p.$originalResidueList[0]$startNumber$mutationResidueList[0]" if($originalProteinVariationLength == 1 && $mutationProteinVariationLength == 1);
	if($originalProteinVariationLength == 0) {
		if(substr($mutationProteinSequence, $proteinVariationPosition - $mutationProteinVariationLength - 1, $mutationProteinVariationLength) eq substr($mutationProteinSequence, $proteinVariationPosition - 1, $mutationProteinVariationLength)) {
			my ($startNumber, $endNumber) = ($proteinVariationPosition - $mutationProteinVariationLength, $proteinVariationPosition - 1);
			return "p.$mutationResidueList[0]${startNumber}dup" if($mutationProteinVariationLength == 1);
			return "p.$mutationResidueList[0]${startNumber}_$mutationResidueList[-1]${endNumber}dup";
		}
		($startNumber, $endNumber) = ($endNumber, $startNumber);
		my ($startResidue, $endResidue) = map {$aaOneToThreeLetter{substr($originalProteinSequence, $_ - 1, 1)}} ($startNumber, $endNumber);
		return "p.$startResidue${startNumber}_$endResidue${endNumber}ins" . join('', @mutationResidueList);
	}
	if($mutationProteinVariationLength == 0) {
		return "p.$originalResidueList[0]${startNumber}del" if($originalProteinVariationLength == 1);
		return "p.$originalResidueList[0]${startNumber}_$originalResidueList[-1]${endNumber}del";
	} else {
		return "p.$originalResidueList[0]${startNumber}delins" . join('', @mutationResidueList) if($originalProteinVariationLength == 1);
		return "p.$originalResidueList[0]${startNumber}_$originalResidueList[-1]${endNumber}delins" . join('', @mutationResidueList);
	}
}

sub getProteinMutation {
	my ($proteinVariationNomenclature) = @_;
	if($proteinVariationNomenclature eq '') {
		return '';
	} elsif($proteinVariationNomenclature eq 'p.=') {
		return 'silent';
	} elsif($proteinVariationNomenclature eq 'p.?') {
		return 'startcodon';
	} elsif($proteinVariationNomenclature =~ /\*$/) {
		return 'nonsense';
	} elsif($proteinVariationNomenclature =~ /ext\*[0-9]+$/) {
		return 'readthrough';
	} elsif($proteinVariationNomenclature =~ /fs\*[0-9]+$/) {
		return 'frameshift';
	} elsif($proteinVariationNomenclature =~ /dup/) {
		return 'inframe';
	} elsif($proteinVariationNomenclature =~ /ins/) {
		return 'inframe';
	} elsif($proteinVariationNomenclature =~ /del/) {
		return 'inframe';
	} else {
		return 'missense';
	}
}

sub getNucleotideVariationNomenclature {
	my ($nucleotideSequence, $nucleotideVariationPosition, $originalNucleotideVariation, $mutationNucleotideVariation, $tokenHash) = @_;
	my $typeLetter = $tokenHash->{'codingRegions'} ne '' ? 'c' : 'r';
	my ($originalNucleotideVariationLength, $mutationNucleotideVariationLength) = map {length} ($originalNucleotideVariation, $mutationNucleotideVariation);
	my ($startNumber, $endNumber) = map {getNucleotideNumber($_, $tokenHash->{'region'}, $tokenHash)} ($nucleotideVariationPosition, $nucleotideVariationPosition + $originalNucleotideVariationLength - 1);
	return "$typeLetter.$startNumber$originalNucleotideVariation>$mutationNucleotideVariation" if($originalNucleotideVariationLength == 1 && $mutationNucleotideVariationLength == 1);
	if($originalNucleotideVariation eq '' || $mutationNucleotideVariation eq '') {
		my $repeatUnit = getRepeatUnit("$originalNucleotideVariation$mutationNucleotideVariation");
		my $repeatCount = 0;
		if(substr($nucleotideSequence, 0, $nucleotideVariationPosition - 1) =~ /($repeatUnit)+$/) {
			$startNumber = getNucleotideNumber($-[0] + 1, $tokenHash->{'region'}, $tokenHash);
			$repeatCount += ($+[0] - $-[0]) / length($repeatUnit);
		}
		if(substr($nucleotideSequence, $nucleotideVariationPosition + $originalNucleotideVariationLength - 1) =~ /^($repeatUnit)+/) {
			$endNumber = getNucleotideNumber($nucleotideVariationPosition + $originalNucleotideVariationLength - 1 + $+[0], $tokenHash->{'region'}, $tokenHash);
			$repeatCount += ($+[0] - $-[0]) / length($repeatUnit);
		}
		my $originalRepeatCount = $repeatCount + $originalNucleotideVariationLength / length($repeatUnit);
		my $mutationRepeatCount = $repeatCount + $mutationNucleotideVariationLength / length($repeatUnit);
		if($originalRepeatCount == 1 && $mutationRepeatCount == 2) {
			return "$typeLetter.${startNumber}dup$mutationNucleotideVariation" if($mutationNucleotideVariationLength == 1);
			return "$typeLetter.${startNumber}_${endNumber}dup$mutationNucleotideVariation";
		} elsif($originalRepeatCount > 0 && $mutationRepeatCount > 0) {
			return "$typeLetter.$startNumber$repeatUnit\[$mutationRepeatCount]";
		}
	}
	if($originalNucleotideVariationLength == 0) {
		($startNumber, $endNumber) = ($endNumber, $startNumber);
		return "$typeLetter.${startNumber}_${endNumber}ins$mutationNucleotideVariation";
	}
	if($mutationNucleotideVariationLength == 0) {
		return "$typeLetter.${startNumber}del$originalNucleotideVariation" if($originalNucleotideVariationLength == 1);
		return "$typeLetter.${startNumber}_${endNumber}del$originalNucleotideVariation";
	} else {
		return "$typeLetter.${startNumber}del${originalNucleotideVariation}ins$mutationNucleotideVariation" if($originalNucleotideVariationLength == 1);
		return "$typeLetter.${startNumber}_${endNumber}del${originalNucleotideVariation}ins$mutationNucleotideVariation";
	}
}

sub getNucleotideNumber {
	my ($position, $region, $tokenHash) = @_;
	if($region eq 'exon') {
		if($tokenHash->{'codingRegions'} ne '') {
			return $position - $tokenHash->{'codingStartInTranscript'} if($position < $tokenHash->{'codingStartInTranscript'});
			return '*' . ($position - $tokenHash->{'codingEndInTranscript'}) if($position > $tokenHash->{'codingEndInTranscript'});
			return $position - $tokenHash->{'codingStartInTranscript'} + 1;
		} else {
			return $position;
		}
	}
	if($region eq 'intron') {
		if($position < $tokenHash->{'intronSequenceLength'} - $position + 1) {
			return sprintf('%s+%d', getNucleotideNumber($tokenHash->{'startInTranscript'}, 'exon', $tokenHash), $position);
		} else {
			return sprintf('%s-%d', getNucleotideNumber($tokenHash->{'endInTranscript'}, 'exon', $tokenHash), $tokenHash->{'intronSequenceLength'} - $position + 1);
		}
	}
}

sub getRepeatUnit {
	my ($sequence) = @_;
	for(my $length = 1; $length <= length($sequence); $length++) {
		my $repeatUnit = substr($sequence, 0, $length);
		return $repeatUnit if($sequence =~ /^($repeatUnit)*$/);
	}
}

sub leftalignIndel {
	my ($sequence, $position, $original, $mutation) = @_;
	if($original eq '' || $mutation eq '') {
		my $indel = "$original$mutation";
		$sequence = substr($sequence, 0, $position - 1);
		$indel = substr($indel, -1, 1, '') . $indel while(substr($indel, -1, 1) eq substr($sequence, -1, 1, '') && $position--);
		$original = $indel if($original ne '');
		$mutation = $indel if($mutation ne '');
	}
	return ($position, $original, $mutation);
}

sub rightalignIndel {
	my ($sequence, $position, $original, $mutation) = @_;
	if($original eq '' || $mutation eq '') {
		my $indel = "$original$mutation";
		$sequence = substr($sequence, $position - 1 + length($original));
		$indel = $indel . substr($indel, 0, 1, '') while(substr($indel, 0, 1) eq substr($sequence, 0, 1, '') && $position++);
		$original = $indel if($original ne '');
		$mutation = $indel if($mutation ne '');
	}
	return ($position, $original, $mutation);
}

sub getExtendedStartEnd {
	my ($chromosome, $position, $refBase, $altBase) = @_;
	my ($start, $end) = ($position, $position + length($refBase) - 1);
	if($refBase eq '' || $altBase eq '') {
		for(my ($indel, $extended) = ("$refBase$altBase", ''); "$indel$extended" =~ /^$extended/;) {
			$end += 1;
			last if($end > $db->length($chromosome));
			$extended .= uc($db->seq($chromosome, $end, $end));
		}
		$end -= 1;
		for(my ($indel, $extended) = ("$refBase$altBase", ''); "$extended$indel" =~ /$extended$/;) {
			$start -= 1;
			last if($start < 1);
			$extended .= uc($db->seq($chromosome, $start, $start));
		}
		$start += 1;
		($start, $end) = ($end, $start) if($end < $start);
	}
	return ($start, $end);
}

sub stripIdentical {
	my ($position, @sequenceList) = @_;
	while(my @baseList = map {substr($_, -1, 1)} @sequenceList) {
		last if(grep {$baseList[$_ - 1] ne $baseList[$_]} 1 .. $#baseList);
		substr($_, -1, 1, '') foreach(@sequenceList);
	}
	while(my @baseList = map {substr($_, 0, 1)} @sequenceList) {
		last if(grep {$baseList[$_ - 1] ne $baseList[$_]} 1 .. $#baseList);
		substr($_, 0, 1, '') foreach(@sequenceList);
		$position += 1;
	}
	return ($position, @sequenceList);
}

sub getChromosomeList {
	my @chromosomeList = ();
	if(my $faiFile = `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null`) {
		chomp($faiFile);
		open(my $reader, $faiFile);
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line);
			push(@chromosomeList, $tokenList[0]);
		}
		close($reader);
	} else {
		open(my $reader, $referenceFastaFile);
		while(my $line = <$reader>) {
			chomp($line);
			push(@chromosomeList, $1) if($line =~ /^>(\S*)/);
		}
		close($reader);
	}
	return @chromosomeList;
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
