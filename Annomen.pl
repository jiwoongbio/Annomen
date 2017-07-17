# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum min max);
use Bio::DB::Fasta;
use Getopt::Long;

GetOptions('c=s' => \(my $codons = ''), 's=s' => \(my $startCodons = 'ATG,CTG,TTG'), 'r' => \(my $checkReferenceVariation = ''));
my ($inputFile, $referenceFastaFile, $annotationTableFile, $transcriptFastaFile, $proteinFastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
my %aaOneToThreeLetter = ('A' => 'Ala', 'B' => 'Asx', 'C' => 'Cys', 'D' => 'Asp', 'E' => 'Glu', 'F' => 'Phe', 'G' => 'Gly', 'H' => 'His', 'I' => 'Ile', 'K' => 'Lys', 'L' => 'Leu', 'M' => 'Met', 'N' => 'Asn', 'P' => 'Pro', 'Q' => 'Gln', 'R' => 'Arg', 'S' => 'Ser', 'T' => 'Thr', 'V' => 'Val', 'W' => 'Trp', 'X' => 'Xxx', 'Y' => 'Tyr', 'Z' => 'Glx');
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
$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} split(/,/, $codons));
sub translate {
	my ($sequence) = @_;
	return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. (length($sequence) / 3) - 1);
}
my @startCodonList = split(/,/, $startCodons);
my %transcriptSequenceHash = ();
{
	my $transcriptId = '';
	open(my $reader, ($transcriptFastaFile =~ /\.gz$/) ? "gzip -dc $transcriptFastaFile |" : $transcriptFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ s/^>// && ($transcriptId = $line));
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
		next if($line =~ s/^>// && ($proteinId = $line));
		$proteinSequenceHash{$proteinId} .= $line;
	}
	close($reader);
	s/U/*/g foreach(values %proteinSequenceHash);
	s/\*?$/*/ foreach(values %proteinSequenceHash);
}
{
	chomp(my @chromosomeList = `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null | xargs cat | cut -f1`);
	chomp(@chromosomeList = `grep '^>' $referenceFastaFile | sed 's/^>//' | sed 's/ .*\$//'`) unless(@chromosomeList);
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
my ($title, @columnList) = ('Annomen', 'mutation', 'region', 'strand', 'spliceDistance', 'geneName', 'transcriptId', 'proteinId', 'nucleotideVariationNomenclature', 'proteinVariationNomenclature');
my $vcf = 0;
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^#/) {
		$vcf = 1 if($line =~ /^##fileformat=VCF/);
		print "$line\n";
	} else {
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $position, $refBase, $altBase) = @tokenList[$vcf ? (0, 1, 3, 4) : (0, 1, 2, 3)];
		my @tokensList = ();
		if(my @altBaseList = grep {/^[A-Za-z]*$/} split(/,/, $altBase)) {
			my @startEndList = map {[getExtendedStartEnd($chromosome, stripIdentical($position, $refBase, $_))]} @altBaseList;
			my ($start, $end) = (min(map {$_->[0]} @startEndList), max(map {$_->[1]} @startEndList));
			my @transcriptIdList = ();
			my %transcriptIdTokenHashListHash = ();
			foreach my $tokenList (getTokenListList($chromosome, $position, $end)) {
				my %tokenHash = ();
				@tokenHash{(getColumnList)} = @$tokenList;
				next if($tokenHash{'end'} < $start);
				if(abs($tokenHash{'spliceDistance'} = getSpliceDistance($start, $end, \%tokenHash)) > 0) {
					setMutationNomenclatures($chromosome, $position, $refBase, \@altBaseList, \%tokenHash);
				}
				$tokenHash{'region'} = getRegion($tokenHash{'strand'}, $start, $end, @tokenHash{'region', 'codingStart', 'codingEnd'});
				my $transcriptId = $tokenHash{'transcriptId'};
				push(@transcriptIdList, $transcriptId) unless(defined($transcriptIdTokenHashListHash{$transcriptId}));
				push(@{$transcriptIdTokenHashListHash{$transcriptId}}, \%tokenHash);
			}
			foreach my $transcriptId (@transcriptIdList) {
				my %tokenHash = ();
				%tokenHash = %{$transcriptIdTokenHashListHash{$transcriptId}->[0]};
				if(scalar(my @tokenHashList = @{$transcriptIdTokenHashListHash{$transcriptId}}) > 1) {
					my @regionList = map {split(/\^/, $_->{'region'})} sort {$a->{'startInTranscript'} <=> $b->{'startInTranscript'}} @tokenHashList;
					$tokenHash{'region'} = join('^', map {$regionList[$_]} grep {$_ == 0 || $regionList[$_ - 1] ne $regionList[$_]} 0 .. $#regionList);
					$tokenHash{'spliceDistance'} = min(map {$_->{'spliceDistance'}} @tokenHashList);
				}
				$tokenHash{'mutation'} = 'junction' if($tokenHash{'spliceDistance'} == 0);
				push(@tokensList, join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}));
			}
		}
		if($vcf) {
			if(@tokensList) {
				my @infoList = grep {$_ ne '.' && $_ !~ /^$title\./} split(/;/, $tokenList[7]);
				for(my $number = 1; $number <= scalar(@tokensList); $number++) {
					my %tokenHash = ();
					@tokenHash{@columnList} = split(/\t/, $tokensList[$number - 1], scalar(@columnList));
					if(scalar(@tokensList) > 1) {
						push(@infoList, map {"$title.${_}_$number=$tokenHash{$_}"} grep {$tokenHash{$_} ne ''} @columnList);
					} else {
						push(@infoList, map {"$title.$_=$tokenHash{$_}"} grep {$tokenHash{$_} ne ''} @columnList);
					}
				}
				$tokenList[7] = join(';', @infoList);
			}
			print join("\t", @tokenList), "\n";
		} else {
			push(@tokensList, join("\t", ('') x scalar(@columnList))) if(scalar(@tokensList) == 0);
			print join("\t", @tokenList, $_), "\n" foreach(@tokensList);
		}
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
	my ($chromosome, $position, $refBase, $altBaseList, $tokenHash) = @_;
	($position, $refBase, my @altBaseList) = stripIdentical($position, $refBase, @$altBaseList);
	my @codingRegionList = map {[split(/\.\./, $_)]} split(/,/, $tokenHash->{'codingRegions'});
	@$tokenHash{'codingStartInTranscript', 'codingEndInTranscript'} = map {$codingRegionList[$_]->[$_]} (0, -1) if(@codingRegionList);
	if($tokenHash->{'region'} eq 'exon') {
		my @mismatchList = map {[split(':', $_, 3)]} split(/,/, $tokenHash->{'mismatch'});
		foreach my $mismatch (@mismatchList) {
			if($mismatch->[0] < $position && $position <= $mismatch->[0] + length($mismatch->[1])) {
				my $head = substr($mismatch->[1], 0, $position - $mismatch->[0]);
				($refBase, @altBaseList) = map {"$head$_"} ($refBase, @altBaseList);
				$position = $mismatch->[0];
			}
			if($mismatch->[0] <= $position + length($refBase) && $position + length($refBase) < $mismatch->[0] + length($mismatch->[1])) {
				my $tail = substr($mismatch->[1], $position + length($refBase) - ($mismatch->[0] + length($mismatch->[1])));
				($refBase, @altBaseList) = map {"$_$tail"} ($refBase, @altBaseList);
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
			($rnaBase, $refBase, @altBaseList) = map {getReverseComplementarySequence($_)} ($rnaBase, $refBase, @altBaseList);
		}
		my $transcriptSequence = $transcriptSequenceHash{my $transcriptId = $tokenHash->{'transcriptId'}};
		my @transcriptProteinVariationNomenclaturesList = ();
		push(@transcriptProteinVariationNomenclaturesList, [getTranscriptProteinVariationNomenclatures($transcriptSequence, $position, $rnaBase, $refBase, $tokenHash, @codingRegionList)]) if($checkReferenceVariation && $rnaBase ne $refBase);
		push(@transcriptProteinVariationNomenclaturesList, [getTranscriptProteinVariationNomenclatures($transcriptSequence, $position, $rnaBase,       $_, $tokenHash, @codingRegionList)]) foreach(@altBaseList);
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
				($intronSequence, $refBase, @altBaseList) = map {getReverseComplementarySequence($_)} ($intronSequence, $refBase, @altBaseList);
			}
			$tokenHash->{'intronSequenceLength'} = length($intronSequence);
#			my @intronVariationNomenclatureList = ();
#			if($tokenHash->{'spliceDistance'} > 0) {
#				push(@intronVariationNomenclatureList, getNucleotideVariationNomenclature($intronSequence,  leftalignIndel($intronSequence, stripIdentical($position, $refBase, $_)), $tokenHash)) foreach(@altBaseList);
#			} else {
#				push(@intronVariationNomenclatureList, getNucleotideVariationNomenclature($intronSequence, rightalignIndel($intronSequence, stripIdentical($position, $refBase, $_)), $tokenHash)) foreach(@altBaseList);
#			}
			my @intronVariationNomenclatureList = map {getNucleotideVariationNomenclature($intronSequence, rightalignIndel($intronSequence, stripIdentical($position, $refBase, $_)), $tokenHash)} @altBaseList;
			$tokenHash->{'nucleotideVariationNomenclature'} = join(',', @intronVariationNomenclatureList);
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
			if((my $index = $tokenHash->{'codonStart'} - 1) > 0) {
				$codingVariationPosition = $codingVariationPosition - $index;
				$originalCodingSequence = substr($originalCodingSequence, $index);
				$mutationCodingSequence = substr($mutationCodingSequence, $index);
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
				if(grep {$_ eq $mutationStartCodon} ($originalStartCodon, @startCodonList)) {
					$proteinVariationNomenclature = getProteinVariationNomenclature($originalProteinSequence, $mutationProteinSequence, $proteinVariationPosition, $originalProteinVariationLength, $mutationProteinVariationLength, $frameshift);
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
		my $partialSequence = substr($sequence, 0, $position - 1);
		$indel = substr($indel, -1, 1, '') . $indel while(substr($indel, -1, 1) eq substr($partialSequence, -1, 1, '') && $position--);
		$original = $indel if($original ne '');
		$mutation = $indel if($mutation ne '');
	}
	return ($position, $original, $mutation);
}

sub rightalignIndel {
	my ($sequence, $position, $original, $mutation) = @_;
	if($original eq '' || $mutation eq '') {
		my $indel = "$original$mutation";
		my $partialSequence = substr($sequence, $position - 1 + length($original));
		$indel = $indel . substr($indel, 0, 1, '') while(substr($indel, 0, 1) eq substr($partialSequence, 0, 1, '') && $position++);
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
			$extended .= uc($db->seq($chromosome, $end += 1, $end));
		}
		$end -= 1;
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

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}
