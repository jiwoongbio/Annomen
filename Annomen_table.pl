# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

exit 1 if(map {(`which $_`) ? () : $_} ('needle', 'stretcher'));

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'c=s' => \(my $chromosomeToChromosomeFile = ''),
	'C=s' => \@codonList,
	'T' => \(my $isTranscriptGenome = ''),
	'i' => \(my $minimumIdentity = 0.8),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl Annomen_table.pl [options] gene.gff genome.fasta transcript.fasta protein.fasta [transcript.gb/] > Annomen_table.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -c FILE  chromosome-to-chromosome file
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 1 (standard)]
         -T       is transcript genome
         -i FLOAT minimum identity between genome and transcript sequences

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
my ($gffFile, $referenceFastaFile, $transcriptFastaFile, $proteinFastaFile, $gbDirectory) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
my @chromosomeList = getChromosomeList();
my %chromosomeIndexHash = map {$chromosomeList[$_] => $_} 0 .. $#chromosomeList;
if($chromosomeToChromosomeFile ne '') {
	open(my $reader, ($chromosomeToChromosomeFile =~ /\.gz$/ ? "gzip -dc $chromosomeToChromosomeFile |" : $chromosomeToChromosomeFile));
	while(my $line = <$reader>) {
		chomp($line);
		my @chromosomeList = split(/\t/, $line, -1);
		if(scalar(my ($index) = grep {defined} @chromosomeIndexHash{@chromosomeList}) == 1) {
			foreach my $chromosome (@chromosomeList) {
				$chromosomeIndexHash{$chromosome} = $index unless(defined($chromosomeIndexHash{$chromosome}));
			}
		}
	}
	close($reader);
}
my %transcriptSequenceHash = ();
if(-r $transcriptFastaFile) {
	my $transcriptId = '';
	open(my $reader, ($transcriptFastaFile =~ /\.gz$/ ? "gzip -dc $transcriptFastaFile |" : $transcriptFastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($transcriptId = $1));
		$transcriptSequenceHash{$transcriptId} .= $line;
	}
	close($reader);
}
my %proteinSequenceHash = ();
if(-r $proteinFastaFile) {
	my $proteinId = '';
	open(my $reader, ($proteinFastaFile =~ /\.gz$/ ? "gzip -dc $proteinFastaFile |" : $proteinFastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($proteinId = $1));
		$proteinSequenceHash{$proteinId} .= $line;
	}
	close($reader);
	s/U/*/g foreach(values %proteinSequenceHash);
	s/\*?$/*/ foreach(values %proteinSequenceHash);
}
my @columnList = ('chromosome', 'start', 'end', 'transcriptId', 'geneName', 'strand', 'region', 'number', 'startInTranscript', 'endInTranscript', 'mismatch', 'proteinId', 'codingStart', 'codingEnd', 'codingRegions', 'frame');
print join("\t", @columnList), "\n";
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1n -k2,2n -k3,3n");
{
	open(my $reader, ($gffFile =~ /\.gz$/ ? "gzip -dc $gffFile | grep -v '^#' | LC_ALL=C sort -t '\t' -k1,1 -k4,4n -k5,5nr |" : "grep -v '^#' $gffFile | LC_ALL=C sort -t '\t' -k1,1 -k4,4n -k5,5nr |"));
	my @idList = ();
	my %tokenHashHash = ();
	my %exonStartEndListHash = ();
	my %proteinIdListHash = ();
	my %codingStartEndListHash = ();
	my %frameHash = ();
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my %tokenHash = ();
		@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line, -1);
		my %attributeHash = ();
		$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^;= ]+)=([^;]+)(;|$)/g);
		$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^;" ]+) +"([^;"]+)"(;|$)/g);
		$tokenHash{'attribute'} = \%attributeHash;
		my @remainIdList = ();
		foreach my $id (@idList) {
			if($tokenHashHash{$id}->{'chromosome'} ne $tokenHash{'chromosome'} || (defined($tokenHashHash{$id}->{'end'}) && $tokenHashHash{$id}->{'end'} < $tokenHash{'start'})) {
				my ($transcriptId, $geneName) = (getAttributeValue($tokenHashHash{$id}, 'transcript_id', 'locus_tag', 'gene_name', 'gene'), getAttributeValue($tokenHashHash{$id}, 'gene_name', 'gene'));
				if(defined(my $proteinIdList = $proteinIdListHash{$id})) {
					foreach my $proteinId (@$proteinIdList) {
						if(scalar(@$proteinIdList) > 1) {
							printTable($tokenHashHash{$id}, "$transcriptId/$proteinId", $geneName, $exonStartEndListHash{$id}, $proteinId, $codingStartEndListHash{$id}->{$proteinId}, $frameHash{$id}->{$proteinId});
						} else {
							printTable($tokenHashHash{$id}, $transcriptId, $geneName, $exonStartEndListHash{$id}, $proteinId, $codingStartEndListHash{$id}->{$proteinId}, $frameHash{$id}->{$proteinId});
						}
					}
				} else {
					printTable($tokenHashHash{$id}, $transcriptId, $geneName, $exonStartEndListHash{$id});
				}
				delete $tokenHashHash{$id};
				delete $exonStartEndListHash{$id};
				delete $proteinIdListHash{$id};
				delete $codingStartEndListHash{$id};
				delete $frameHash{$id};
			} else {
				push(@remainIdList, $id);
			}
		}
		@idList = @remainIdList if(scalar(@remainIdList) < scalar(@idList));
		if($tokenHash{'feature'} eq 'exon' || $tokenHash{'feature'} eq 'CDS') {
			my $id = $tokenHash{'attribute'}->{'Parent'};
			unless(defined($id)) {
				$id = getAttributeValue(\%tokenHash, 'transcript_id', 'locus_tag');
				unless(defined($tokenHashHash{$id})) {
					$tokenHashHash{$id}->{'attribute'}->{'transcript_id'} = $id;
					@{$tokenHashHash{$id}}{'chromosome', 'strand'} = @tokenHash{'chromosome', 'strand'};
					push(@idList, $id);
				}
				unless(defined($tokenHashHash{$id}->{'attribute'}->{'gene_name'})) {
					if(defined(my $geneName = getAttributeValue(\%tokenHash, 'gene_name', 'gene'))) {
						$tokenHashHash{$id}->{'attribute'}->{'gene_name'} = $geneName;
					}
				}
			}
			if($tokenHash{'feature'} eq 'exon') {
				push(@{$exonStartEndListHash{$id}}, [@tokenHash{'start', 'end'}]);
			}
			if($tokenHash{'feature'} eq 'CDS') {
				my $proteinId = getAttributeValue(\%tokenHash, 'protein_id');
				unless(defined($proteinId)) {
					($proteinId) = grep {defined($proteinSequenceHash{$_})} getAttributeValue(\%tokenHash, 'transcript_id', 'locus_tag', 'gene_name', 'gene');
				}
				$proteinId = $id unless(defined($proteinId));
				push(@{$proteinIdListHash{$id}}, $proteinId) unless(defined($codingStartEndListHash{$id}->{$proteinId}));
				push(@{$codingStartEndListHash{$id}->{$proteinId}}, [@tokenHash{'start', 'end'}]);
				if($tokenHash{'strand'} eq '+') {
					$frameHash{$id}->{$proteinId} = $tokenHash{'frame'} unless(defined($frameHash{$id}->{$proteinId}));
				}
				if($tokenHash{'strand'} eq '-') {
					$frameHash{$id}->{$proteinId} = $tokenHash{'frame'};
				}
			}
		} elsif(defined(my $id = $tokenHash{'attribute'}->{'ID'})) {
			$tokenHashHash{$id} = \%tokenHash;
			push(@idList, $id);
		}
	}
	foreach my $id (@idList) {
		my ($transcriptId, $geneName) = (getAttributeValue($tokenHashHash{$id}, 'transcript_id', 'locus_tag'), getAttributeValue($tokenHashHash{$id}, 'gene_name', 'gene'));
		if(defined(my $proteinIdList = $proteinIdListHash{$id})) {
			foreach my $proteinId (@$proteinIdList) {
				if(scalar(@$proteinIdList) > 1) {
					printTable($tokenHashHash{$id}, "$transcriptId/$proteinId", $geneName, $exonStartEndListHash{$id}, $proteinId, $codingStartEndListHash{$id}->{$proteinId}, $frameHash{$id}->{$proteinId});
				} else {
					printTable($tokenHashHash{$id}, $transcriptId, $geneName, $exonStartEndListHash{$id}, $proteinId, $codingStartEndListHash{$id}->{$proteinId}, $frameHash{$id}->{$proteinId});
				}
			}
		} else {
			printTable($tokenHashHash{$id}, $transcriptId, $geneName, $exonStartEndListHash{$id});
		}
		delete $tokenHashHash{$id};
		delete $exonStartEndListHash{$id};
		delete $proteinIdListHash{$id};
		delete $codingStartEndListHash{$id};
		delete $frameHash{$id};
	}
	@idList = ();
	close($reader);

	sub getAttributeValue {
		my ($tokenHash, @attributeList) = @_;
		foreach my $attribute (@attributeList) {
			if(defined(my $value = $tokenHash->{'attribute'}->{$attribute})) {
				return $value;
			}
		}
		foreach my $attribute (@attributeList) {
			if(defined(my $parent = $tokenHash->{'attribute'}->{'Parent'})) {
				return getAttributeValue($tokenHashHash{$parent}, @attributeList);
			}
		}
		return undef;
	}
}
close($writer);
{
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		$tokenHash{'chromosome'} = $chromosomeList[$tokenHash{'chromosome'}];
		print join("\t", @tokenHash{@columnList}), "\n";
	}
}
close($reader);
waitpid($pid, 0);
unless(-r $transcriptFastaFile) {
	open(my $writer, ($transcriptFastaFile =~ /\.gz$/ ? "| gzip > $transcriptFastaFile" : "> $transcriptFastaFile"));
	foreach my $transcriptId (sort keys %transcriptSequenceHash) {
		my $transcriptSequence = $transcriptSequenceHash{$transcriptId};
		print $writer ">$transcriptId\n";
		for(my $index = 0; $index < length($transcriptSequence); $index += 80) {
			print $writer substr($transcriptSequence, $index, 80), "\n";
		}
	}
	close($writer);
}
unless(-r $proteinFastaFile) {
	open(my $writer, ($proteinFastaFile =~ /\.gz$/ ? "| gzip > $proteinFastaFile" : "> $proteinFastaFile"));
	foreach my $proteinId (sort keys %proteinSequenceHash) {
		my $proteinSequence = $proteinSequenceHash{$proteinId};
		print $writer ">$proteinId\n";
		for(my $index = 0; $index < length($proteinSequence); $index += 80) {
			print $writer substr($proteinSequence, $index, 80), "\n";
		}
	}
	close($writer);
}

sub printTable {
	my ($tokenHash, $transcriptId, $geneName, $exonStartEndList, $proteinId, $codingStartEndList, $frame) = @_;
	return unless(defined($exonStartEndList) || defined($codingStartEndList));
	unless(defined($exonStartEndList)) {
		$exonStartEndList = [map {[@$_]} @$codingStartEndList];
	}
	if($isTranscriptGenome) {
		$exonStartEndList->[0]->[0] = 1;
		$exonStartEndList->[-1]->[1] = $db->length($tokenHash->{'chromosome'});
	}
	my ($chromosome, $strand) = @$tokenHash{'chromosome', 'strand'};
	return unless(defined(my $chromosomeIndex = $chromosomeIndexHash{$chromosome}));
	return unless(defined($transcriptId));
	my $exonCount = scalar(my @exonStartEndList = @$exonStartEndList);
	my $exonSequence = '';
	my @exonPositionList = ();
	foreach(@exonStartEndList) {
		my ($exonStart, $exonEnd) = @$_;
		$exonSequence .= $db->seq($chromosomeList[$chromosomeIndex], $exonStart, $exonEnd);
		push(@exonPositionList, $exonStart .. $exonEnd);
	}
	$exonSequence =~ tr/a-z/A-Z/;
	if($strand eq '-') {
		($exonSequence = reverse($exonSequence)) =~ tr/ACGT/TGCA/;
		@exonPositionList = reverse(@exonPositionList);
	}
	my @seq_objectList = ();
	if(defined($gbDirectory)) {
		if(-e "$gbDirectory/$transcriptId.gb") {
			my $seqInput_object = Bio::SeqIO->new(-file => "$gbDirectory/$transcriptId.gb", -format => 'genbank');
			while(my $seq_object = $seqInput_object->next_seq()) {
				push(@seq_objectList, $seq_object);
			}
		} elsif($transcriptId =~ /^(.+)\.[0-9]+$/ && -e "$gbDirectory/$1.gb") {
			my $seqInput_object = Bio::SeqIO->new(-file => "$gbDirectory/$1.gb", -format => 'genbank');
			while(my $seq_object = $seqInput_object->next_seq()) {
				push(@seq_objectList, $seq_object);
			}
		} else {
			print STDERR join("\t", $transcriptId, 'no gb'), "\n";
		}
	} else {
		push(@seq_objectList, '');
	}
	foreach my $seq_object (@seq_objectList) {
		my $transcriptSequence = $transcriptSequenceHash{$transcriptId};
		unless(defined($transcriptSequence)) {
			if(-r $transcriptFastaFile) {
				print STDERR join("\t", $transcriptId, 'no transcript sequence'), "\n";
				next;
			} else {
				$transcriptSequence = $transcriptSequenceHash{$transcriptId} = $exonSequence;
			}
		}
		if($seq_object) {
			$transcriptId = join('.', $seq_object->display_id, $seq_object->version);
			if($transcriptSequence ne $seq_object->seq) {
				print STDERR join("\t", $transcriptId, 'different transcript sequences', $transcriptSequence, $seq_object->seq), "\n";
				next;
			}
		}
		my ($exonAlignment, $transcriptAlignment, $identity) = ($exonSequence, $transcriptSequence, 1);
		if($exonAlignment ne $transcriptAlignment) {
			if($transcriptAlignment =~ /$exonAlignment/) {
				$identity = length($exonAlignment) / length($transcriptAlignment);
				$exonAlignment = ('-' x index($transcriptAlignment, $exonAlignment)) . $exonAlignment;
				$exonAlignment .= '-' x length($1) if(substr($transcriptAlignment, length($exonAlignment)) =~ /^(.*[^A])A*$/);
			} elsif($exonAlignment =~ /$transcriptAlignment/) {
				$identity = length($transcriptAlignment) / length($exonAlignment);
				$transcriptAlignment = ('-' x index($exonAlignment, $transcriptAlignment)) . $transcriptAlignment;
				$transcriptAlignment .= '-' x (length($exonAlignment) - length($transcriptAlignment));
			} else {
				($exonAlignment, $transcriptAlignment, $identity) = needle($exonSequence, $transcriptSequence, 'EDNAFULL', 20, 0);
				($exonAlignment, $transcriptAlignment, $identity) = stretcher($exonSequence, $transcriptSequence, 'EDNAFULL', 20, 0) if($identity eq '');
				$exonAlignment =~ s/-+$//;
				if(length($exonAlignment) < length($transcriptAlignment)) {
					$exonAlignment .= '-' x length($1) if(substr($transcriptAlignment, length($exonAlignment)) =~ /^(.*[^A])A*$/);
					$transcriptAlignment = substr($transcriptAlignment, 0, length($exonAlignment));
				} else {
					$transcriptAlignment .= '-' x (length($exonAlignment) - length($transcriptAlignment));
				}
				if($strand eq '+') {
					$exonAlignment = leftalignGaps($exonAlignment, $transcriptAlignment);
					$transcriptAlignment = leftalignGaps($transcriptAlignment, $exonAlignment);
				}
				if($strand eq '-') {
					$exonAlignment = rightalignGaps($exonAlignment, $transcriptAlignment);
					$transcriptAlignment = rightalignGaps($transcriptAlignment, $exonAlignment);
				}
			}
		}
		print STDERR join("\t", $transcriptId, 'identity', $identity), "\n";
		next if($identity < $minimumIdentity);
		my %exon2transcriptPositionHash = ();
		my %transcript2exonPositionHash = ();
		my @mismatchList = ();
		for(my ($index, $exonIndex, $transcriptIndex) = (0, 0, 0); $index < length($exonAlignment) && $index < length($transcriptAlignment); $index++) {
			my $exonBase = substr($exonAlignment, $index, 1);
			my $transcriptBase = substr($transcriptAlignment, $index, 1);
			if($exonBase ne '-' && $transcriptBase ne '-') {
				$exon2transcriptPositionHash{$exonPositionList[$exonIndex]} = $transcriptIndex + 1;
				$transcript2exonPositionHash{$transcriptIndex + 1} = $exonPositionList[$exonIndex];
			}
			if($exonBase ne $transcriptBase) {
				my $exonPosition = $exonPositionList[$exonIndex - ($strand eq '-' && $exonBase eq '-' ? 1 : 0)];
				push(@mismatchList, [$exonPosition, $exonBase, $transcriptBase, $index]) if(defined($exonPosition));
			}
			$exonIndex++ if($exonBase ne '-');
			$transcriptIndex++ if($transcriptBase ne '-');
		}
		if(@mismatchList) {
			$_->[1] =~ s/-//g foreach(@mismatchList);
			$_->[2] =~ s/-//g foreach(@mismatchList);
			if($strand eq '-') {
				$_->[1] =~ tr/ACGT/TGCA/ foreach(@mismatchList);
				$_->[2] =~ tr/ACGT/TGCA/ foreach(@mismatchList);
				@mismatchList = reverse(@mismatchList);
			}
		}
		my @cdsList = ();
		if($seq_object) {
			foreach my $feat_object (grep {$_->primary_tag eq 'CDS'} $seq_object->get_SeqFeatures) {
				my ($proteinId) = $feat_object->get_tag_values('protein_id');
				my ($codonStart) = $feat_object->get_tag_values('codon_start');
				my $codingRegions = join(',', map {join('..', $_->start, $_->end)} $feat_object->location->each_Location);
				if($codingRegions =~ s/([0-9]+)$//) {
					if(translate(substr($transcriptSequence, $1 - 3, 3)) ne '*' && translate(substr($transcriptSequence, $1, 3)) eq '*') {
						$codingRegions .= $1 + 3;
					} else {
						$codingRegions .= $1;
					}
				}
				my @codingTranscriptPositionList = eval($codingRegions);
				my @codingExonPositionList = sort {$a <=> $b} grep {defined} @transcript2exonPositionHash{@codingTranscriptPositionList};
				my $frame = $codonStart - 1;
				push(@cdsList, [$proteinId, @codingExonPositionList[0, -1], $codingRegions, $frame]);
			}
		} else {
			if(defined($proteinId)) {
				my @codingStartEndTranscriptPositionList = @exon2transcriptPositionHash{$codingStartEndList->[0]->[0], $codingStartEndList->[-1]->[1]};
				if(defined($codingStartEndTranscriptPositionList[0]) && defined($codingStartEndTranscriptPositionList[1])) {
					my $codingRegions = join('..', sort {$a <=> $b} @codingStartEndTranscriptPositionList);
					if($codingRegions =~ s/([0-9]+)$//) {
						if(translate(substr($transcriptSequence, $1 - 3, 3)) ne '*' && translate(substr($transcriptSequence, $1, 3)) eq '*') {
							$codingRegions .= $1 + 3;
						} else {
							$codingRegions .= $1;
						}
					}
					my @codingTranscriptPositionList = eval($codingRegions);
					my @codingExonPositionList = sort {$a <=> $b} grep {defined} @transcript2exonPositionHash{@codingTranscriptPositionList};
					push(@cdsList, [$proteinId, @codingExonPositionList[0, -1], $codingRegions, $frame]);
				} else {
					print STDERR join("\t", $transcriptId, $proteinId, 'unmatched coding start/end position'), "\n";
				}
			}
		}
		my %cdsHash = ();
		foreach my $cds (@cdsList) {
			my ($proteinId, $codingStart, $codingEnd, $codingRegions, $frame) = @$cds;
			my $codingSequence = join('', map {substr($transcriptSequence, $_->[0] - 1, $_->[-1] - ($_->[0] - 1))} map {[split(/\.\./, $_, 2)]} split(/,/, $codingRegions));
			my $translateSequence = translate(substr($codingSequence, $frame));
			my $proteinSequence = $proteinSequenceHash{$proteinId};
			unless(defined($proteinSequence)) {
				if(-r $proteinFastaFile) {
					print STDERR join("\t", $transcriptId, $proteinId, 'no protein sequence'), "\n";
					$cdsHash{$cds} = 1;
					next;
				} else {
					$proteinSequence = $proteinSequenceHash{$proteinId} = $translateSequence;
				}
			}
			if($frame > 0) {
				(my $translateTerminateSequence = $translateSequence) =~ s/\*.*$//;
				if($proteinSequence =~ /^X$translateTerminateSequence/) {
					$translateSequence = "X$translateSequence";
					$cds->[4] = $frame = $frame - 3;
				}
			}
			if($translateSequence ne $proteinSequence) {
				for(my $index = 0; $index < length($proteinSequence) || $index < length($proteinSequence); $index++) {
					my $proteinAA = $index < length($proteinSequence) ? substr($proteinSequence, $index, 1) : '';
					my $translateAA = $index < length($translateSequence) ? substr($translateSequence, $index, 1) : '';
					print STDERR join("\t", $transcriptId, $proteinId, 'different amino acids', $index + 1, $proteinAA, $translateAA), "\n" if($translateAA ne $proteinAA);
				}
			}
		}
		@cdsList = grep {!defined($cdsHash{$_})} @cdsList;
		for(my $index = 0; $index < $exonCount; $index++) {
			my $number = $index + 1;
			$number = $exonCount - $index if($strand eq '-');
			my ($region, $start, $end) = ('exon', @{$exonStartEndList[$index]});
			next if(!defined($exon2transcriptPositionHash{$start}));
			next if(!defined($exon2transcriptPositionHash{$end}));
			my %tokenHash = ();
			@tokenHash{'chromosome', 'start', 'end', 'transcriptId', 'geneName', 'strand', 'region', 'number'} = ($chromosomeIndex, $start, $end, $transcriptId, $geneName, $strand, $region, $number);
			@tokenHash{'startInTranscript', 'endInTranscript'} = sort {$a <=> $b} @exon2transcriptPositionHash{$start, $end};
			if(my @includedMismatchList = grep {!($_->[0] == $start && $_->[1] eq '')} grep {$start <= $_->[0] && $_->[0] <= $end} @mismatchList) {
				my @indexList = grep {abs($includedMismatchList[$_]->[3] - $includedMismatchList[$_ - 1]->[3]) > 1} 1 .. $#includedMismatchList;
				@indexList = (0, @indexList, scalar(@includedMismatchList));
				@includedMismatchList = map {[$_->[0]->[0], join('', map {$_->[1]} @$_), join('', map {$_->[2]} @$_)]} map {[@includedMismatchList[$indexList[$_ - 1] .. $indexList[$_] - 1]]} 1 .. $#indexList;
				$tokenHash{'mismatch'} = join(',', map {join(':', @$_)} @includedMismatchList);
			}
			if(@cdsList) {
				foreach(@cdsList) {
					@tokenHash{'proteinId', 'codingStart', 'codingEnd', 'codingRegions', 'frame'} = @$_;
					print $writer join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
				}
			} else {
				print $writer join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
			}
		}
		for(my $index = 1; $index < $exonCount; $index++) {
			my $number = $index;
			$number = $exonCount - $index if($strand eq '-');
			my ($region, $start, $end) = ('intron', $exonStartEndList[$index - 1]->[1] + 1, $exonStartEndList[$index]->[0] - 1);
			next if($start > $end);
			next if(!defined($exon2transcriptPositionHash{$start - 1}));
			next if(!defined($exon2transcriptPositionHash{$end + 1}));
			my %tokenHash = ();
			@tokenHash{'chromosome', 'start', 'end', 'transcriptId', 'geneName', 'strand', 'region', 'number'} = ($chromosomeIndex, $start, $end, $transcriptId, $geneName, $strand, $region, $number);
			@tokenHash{'startInTranscript', 'endInTranscript'} = sort {$a <=> $b} @exon2transcriptPositionHash{$start - 1, $end + 1};
			if(@cdsList) {
				foreach(@cdsList) {
					@tokenHash{'proteinId', 'codingStart', 'codingEnd', 'codingRegions', 'frame'} = @$_;
					print $writer join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
				}
			} else {
				print $writer join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}), "\n";
			}
		}
	}
}

sub leftalignGaps {
	my ($alignA, $alignB) = @_;
	my @gapStartEndList = ();
	push(@gapStartEndList, [$-[0] + 1, $+[0] - 1]) while($alignA =~ /[^-]-+[^-]/g);
	foreach(@gapStartEndList) {
		my ($gapStart, $gapEnd) = @$_;
		while(substr($alignA, $gapStart - 1, 1) eq substr($alignB, $gapEnd - 1, 1)) {
			$alignA = join('', map {substr($alignA, $_->[0], $_->[1] - $_->[0])} ([0, $gapStart - 1], [$gapStart, $gapEnd], [$gapStart - 1, $gapStart], [$gapEnd, length($alignA)]));
			($gapStart, $gapEnd) = map {$_ - 1} ($gapStart, $gapEnd);
		}
	}
	return $alignA;
}

sub rightalignGaps {
	my ($alignA, $alignB) = @_;
	my @gapStartEndList = ();
	push(@gapStartEndList, [$-[0] + 1, $+[0] - 1]) while($alignA =~ /[^-]-+[^-]/g);
	foreach(reverse @gapStartEndList) {
		my ($gapStart, $gapEnd) = @$_;
		while(substr($alignA, $gapEnd, 1) eq substr($alignB, $gapStart, 1)) {
			$alignA = join('', map {substr($alignA, $_->[0], $_->[1] - $_->[0])} ([0, $gapStart], [$gapEnd, $gapEnd + 1], [$gapStart, $gapEnd], [$gapEnd + 1, length($alignA)]));
			($gapStart, $gapEnd) = map {$_ + 1} ($gapStart, $gapEnd);
		}
	}
	return $alignA;
}

sub needle {
	my ($seqA, $seqB, $datafile, $gapOpen, $gapExtend) = @_;
	$gapOpen = '' unless(defined($gapOpen));
	$gapExtend = '' unless(defined($gapExtend));
	my $prefix = "$temporaryDirectory/needle.$hostname.$$";
	system("rm -fr $prefix.*");
	{
		open(my $writer, "> $prefix.asequence");
		print $writer ">seqA\n";
		print $writer "$seqA\n";
		close($writer);
	}
	{
		open(my $writer, "> $prefix.bsequence");
		print $writer ">seqB\n";
		print $writer "$seqB\n";
		close($writer);
	}
	print STDERR `echo -en '$gapOpen\\n$gapExtend\\n' | timeout 600 needle -asequence $prefix.asequence -bsequence $prefix.bsequence -outfile $prefix.outfile -datafile $datafile`;
	print STDERR "\n";
	my ($alignA, $alignB, $identity) = ('', '', '');
	if(-e "$prefix.outfile") {
		open(my $reader, "$prefix.outfile");
		while(my $line = <$reader>) {
			chomp($line);
			$identity = $1 / $2 if($line =~ /^# Identity: +([0-9]+)\/([0-9]+) \( *[0-9.]+%\)$/);
			$alignA .= $1 if($line =~ /^seqA +[0-9]+ +([^ ]+) +[0-9]+$/);
			$alignB .= $1 if($line =~ /^seqB +[0-9]+ +([^ ]+) +[0-9]+$/);
		}
		close($reader);
	}
	system("rm -fr $prefix.*");
	return ($alignA, $alignB, $identity);
}

sub stretcher {
	my ($seqA, $seqB, $datafile, $gapOpen, $gapExtend) = @_;
	$gapOpen = defined($gapOpen) && $gapOpen ne '' ? "-gapopen $gapOpen" : '';
	$gapExtend = defined($gapExtend) && $gapExtend ne '' ? "-gapextend $gapExtend" : '';
	my $prefix = "$temporaryDirectory/stretcher.$hostname.$$";
	system("rm -fr $prefix.*");
	{
		open(my $writer, "> $prefix.asequence");
		print $writer ">seqA\n";
		print $writer "$seqA\n";
		close($writer);
	}
	{
		open(my $writer, "> $prefix.bsequence");
		print $writer ">seqB\n";
		print $writer "$seqB\n";
		close($writer);
	}
	print STDERR `timeout 600 stretcher -asequence $prefix.asequence -bsequence $prefix.bsequence -outfile $prefix.outfile -datafile $datafile $gapOpen $gapExtend`;
	print STDERR "\n";
	my ($alignA, $alignB, $identity) = ('', '', '');
	if(-e "$prefix.outfile") {
		open(my $reader, "$prefix.outfile");
		while(my $line = <$reader>) {
			chomp($line);
			$identity = $1 / $2 if($line =~ /^# Identity: +([0-9]+)\/([0-9]+) \( *[0-9.]+%\)$/);
			$alignA .= $1 if($line =~ /^ *seqA (.+)$/);
			$alignB .= $1 if($line =~ /^ *seqB (.+)$/);
		}
		close($reader);
	}
	system("rm -fr $prefix.*");
	return ($alignA, $alignB, $identity);
}

sub getChromosomeList {
	my @chromosomeList = ();
	if(my $faiFile = `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null`) {
		chomp($faiFile);
		open(my $reader, $faiFile);
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line, -1);
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
