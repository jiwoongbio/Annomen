# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

exit 1 if(map {(`which $_`) ? () : $_} ('needle', 'stretcher'));

my @codonList = ();
GetOptions('c=s' => \@codonList);
my ($gtfFile, $referenceFastaFile, $transcriptFastaFile, $proteinFastaFile, $gbDirectory) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
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
$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} map {split(/,/, $_)} @codonList);
sub translate {
	my ($sequence) = @_;
	return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. (length($sequence) / 3) - 1);
}
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
my @columnList = ('chromosome', 'start', 'end', 'transcriptId', 'geneName', 'strand', 'region', 'number', 'startInTranscript', 'endInTranscript', 'mismatch', 'proteinId', 'codingStart', 'codingEnd', 'codingRegions', 'frame');
print join("\t", @columnList), "\n";
my @chromosomeList = getChromosomeList();
foreach my $chromosome (@chromosomeList) {
	my %transcriptGeneNameHash = ();
	my %transcriptStrandHash = ();
	my %transcriptExonStartEndListHash = ();

	my %transcriptProteinHash = ();
	my %transcriptCodingStartHash = ();
	my %transcriptCodingEndHash = ();
	my %transcriptFrameHash = ();

	open(my $reader, "awk -F'\t' '(\$1 == \"$chromosome\")' $gtfFile | sort --field-separator='\t' -k4,4n -k5,5n |");
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my %tokenHash = ();
		@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
		my %attributeHash = ();
		$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);
		if($tokenHash{'feature'} eq 'exon' || $tokenHash{'feature'} eq 'CDS') {
			my $transcriptId = $attributeHash{'transcript_id'};
			$transcriptGeneNameHash{$transcriptId} = $_ if(defined($_ = $attributeHash{'gene_name'}));
			$transcriptStrandHash{$transcriptId} = $tokenHash{'strand'};
			if($tokenHash{'feature'} eq 'exon') {
				push(@{$transcriptExonStartEndListHash{$transcriptId}}, [@tokenHash{'start', 'end'}]);
			}
			$transcriptProteinHash{$transcriptId} = $_ if(defined($_ = $attributeHash{'protein_id'}));
			if($tokenHash{'feature'} eq 'CDS') {
				$transcriptCodingStartHash{$transcriptId} = $tokenHash{'start'} unless(defined($transcriptCodingStartHash{$transcriptId}));
				$transcriptCodingEndHash{$transcriptId} = $tokenHash{'end'};
				if($tokenHash{'strand'} eq '+') {
					$transcriptFrameHash{$transcriptId} = $tokenHash{'frame'} unless(defined($transcriptFrameHash{$transcriptId}));
				}
				if($tokenHash{'strand'} eq '-') {
					$transcriptFrameHash{$transcriptId} = $tokenHash{'frame'};
				}
			}
		}
	}
	close($reader);
	open(my $writer, "| sort --field-separator='\t' -k2,2n -k3,3n");
	foreach my $transcriptId (sort keys %transcriptExonStartEndListHash) {
		printTable($writer, $chromosome, $transcriptId, $transcriptGeneNameHash{$transcriptId}, $transcriptStrandHash{$transcriptId}, $transcriptProteinHash{$transcriptId}, $transcriptCodingStartHash{$transcriptId}, $transcriptCodingEndHash{$transcriptId}, $transcriptFrameHash{$transcriptId}, @{$transcriptExonStartEndListHash{$transcriptId}});
	}
	close($writer);
}

sub printTable {
	my ($writer, $chromosome, $transcriptId, $geneName, $strand, $proteinId, $codingStart, $codingEnd, $frame, @exonStartEndList) = @_;
	$transcriptId =~ s/_dup[0-9]+$//;
	my $exonCount = scalar(@exonStartEndList);
	my $exonSequence = '';
	my @exonPositionList = ();
	foreach(@exonStartEndList) {
		my ($exonStart, $exonEnd) = @$_;
		$exonSequence .= $db->seq($chromosome, $exonStart, $exonEnd);
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
		if($seq_object) {
			$transcriptId = join('.', $seq_object->display_id, $seq_object->version);
			$transcriptSequence = $seq_object->seq;
		}
		unless(defined($transcriptSequenceHash{$transcriptId})) {
			print STDERR join("\t", $transcriptId, 'no transcript sequence'), "\n";
			next;
		}
		if($transcriptSequenceHash{$transcriptId} ne $transcriptSequence) {
			print STDERR join("\t", $transcriptId, 'different transcript sequences', $transcriptSequenceHash{$transcriptId}, $transcriptSequence), "\n";
			next;
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
				($exonAlignment, $transcriptAlignment, $identity) = needle($exonSequence, $transcriptSequence, 20, 0);
				($exonAlignment, $transcriptAlignment, $identity) = stretcher($exonSequence, $transcriptSequence, 20, 0) if($identity == 0);
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
				my ($proteinSequence) = $feat_object->get_tag_values('translation');
				$proteinSequence =~ s/U/*/g;
				$proteinSequence =~ s/\*?$/*/;
				unless(defined($proteinSequenceHash{$proteinId})) {
					print STDERR join("\t", $transcriptId, $proteinId, 'no protein sequence'), "\n";
					next;
				}
				if($proteinSequenceHash{$proteinId} ne $proteinSequence) {
					print STDERR join("\t", $transcriptId, $proteinId, 'different protein sequences', $proteinSequenceHash{$proteinId}, $proteinSequence), "\n";
					next;
				}
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
			if(defined($proteinId) && defined($codingStart) && defined($codingEnd)) {
				unless(defined($proteinSequenceHash{$proteinId})) {
					print STDERR join("\t", $transcriptId, $proteinId, 'no protein sequence'), "\n";
					next;
				}
				my $codingRegions = join('..', sort {$a <=> $b} @exon2transcriptPositionHash{$codingStart, $codingEnd});
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
			}
		}
		foreach my $cds (@cdsList) {
			my ($proteinId, $codingStart, $codingEnd, $codingRegions, $frame) = @$cds;
			my $proteinSequence = $proteinSequenceHash{$proteinId};
			my $codingSequence = join('', map {substr($transcriptSequence, $_->[0] - 1, $_->[-1] - ($_->[0] - 1))} map {[split(/\.\./, $_)]} split(/,/, $codingRegions));
			my $translateSequence = translate(substr($codingSequence, $frame));
			if($frame > 0 && $proteinSequence =~ /^X$translateSequence/) {
				$translateSequence = "X$translateSequence";
				$cds->[4] = $frame = $frame - 3;
			}
			if($translateSequence ne $proteinSequence) {
				for(my $index = 0; $index < length($proteinSequence) || $index < length($proteinSequence); $index++) {
					my $proteinAA = $index < length($proteinSequence) ? substr($proteinSequence, $index, 1) : '';
					my $translateAA = $index < length($translateSequence) ? substr($translateSequence, $index, 1) : '';
					print STDERR join("\t", $transcriptId, $proteinId, 'different amino acids', $index + 1, $translateAA, $proteinAA), "\n" if($translateAA ne $proteinAA);
				}
			}
		}
		for(my $index = 0; $index < $exonCount; $index++) {
			my $number = $index + 1;
			$number = $exonCount - $index if($strand eq '-');
			my ($region, $start, $end) = ('exon', @{$exonStartEndList[$index]});
			next if(!defined($exon2transcriptPositionHash{$start}));
			next if(!defined($exon2transcriptPositionHash{$end}));
			my %tokenHash = ();
			@tokenHash{'chromosome', 'start', 'end', 'transcriptId', 'geneName', 'strand', 'region', 'number'} = ($chromosome, $start, $end, $transcriptId, $geneName, $strand, $region, $number);
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
			@tokenHash{'chromosome', 'start', 'end', 'transcriptId', 'geneName', 'strand', 'region', 'number'} = ($chromosome, $start, $end, $transcriptId, $geneName, $strand, $region, $number);
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
	my ($seqA, $seqB, $gapOpen, $gapExtend) = @_;
	$gapOpen = '' unless(defined($gapOpen));
	$gapExtend = '' unless(defined($gapExtend));
	my $temporaryDirectory = $ENV{'TMPDIR'};
	$temporaryDirectory = '/tmp' unless($temporaryDirectory);
	my $prefix = "$temporaryDirectory/$$.needle";
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
	print STDERR `bash -c "echo -en '$gapOpen\\n$gapExtend\\n'" | timeout 600 needle -asequence $prefix.asequence -bsequence $prefix.bsequence -outfile $prefix.outfile -datafile EDNAFULL`;
	print STDERR "\n";
	my ($alignA, $alignB, $identity) = ('', '', 0);
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
	my ($seqA, $seqB, $gapOpen, $gapExtend) = @_;
	$gapOpen = defined($gapOpen) ? "-gapopen $gapOpen" : '';
	$gapExtend = defined($gapExtend) ? "-gapextend $gapExtend" : '';
	my $temporaryDirectory = $ENV{'TMPDIR'};
	$temporaryDirectory = '/tmp' unless($temporaryDirectory);
	my $prefix = "$temporaryDirectory/$$.stretcher";
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
	print STDERR `timeout 600 stretcher -asequence $prefix.asequence -bsequence $prefix.bsequence -outfile $prefix.outfile -datafile EDNAFULL $gapOpen $gapExtend`;
	print STDERR "\n";
	my ($alignA, $alignB, $identity) = ('', '', 0);
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
