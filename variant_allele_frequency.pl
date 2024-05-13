# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max);
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 1000),
	'c' => \(my $copyBam = ''),
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'S=s' => \(my $stranded = ''),
	'doNotUseSamModule' => \(my $doNotUseSamModule = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl variant_allele_frequency.pl [options] variant.txt reference.fasta sample.bam > sample.variant_allele_frequency.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -c       copy bam file to temporary directory
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -f INT   include flag [$includeFlag]
         -F INT   exclude flag [$excludeFlag]
         -S STR   stranded, "forward", "f", "reverse", or "r"

EOF
}
my $temporaryPrefix = "$temporaryDirectory/$hostname.$$";
system("rm -f $temporaryPrefix.*");
{
	my %pidHash = ();
	my @pidList = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			push(@pidList, $pid);
		} else {
			open($writer, "> $temporaryPrefix.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(@pidList) >= $number) {
			$pidHash{my $pid = wait()} = 1;
			if($pidHash{$pidList[0]}) {
				my $pid = shift(@pidList);
				open(my $reader, "$temporaryPrefix.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryPrefix.$pid");
				delete $pidHash{$pid};
			}
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
my ($variantFile, $referenceFastaFile, $bamFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
my $samModule = $doNotUseSamModule ? '' : eval {require Bio::DB::Sam; 1;};
if($copyBam) {
	system("cp --preserve=timestamps $bamFile $temporaryPrefix.bam");
	if(-r "$bamFile.bai") {
		system("cp --preserve=timestamps $bamFile.bai $temporaryPrefix.bam.bai");
	} else {
		(my $baiFile = $bamFile) =~ s/\.bam/.bai/;
		if(-r $baiFile) {
			system("cp --preserve=timestamps $baiFile $temporaryPrefix.bam.bai");
		} else {
			system("samtools index $temporaryPrefix.bam");
		}
	}
	$bamFile = "$temporaryPrefix.bam";
} else {
	if(-r "$bamFile.bai") {
	} else {
		(my $baiFile = $bamFile) =~ s/\.bam/.bai/;
		if(-r $baiFile) {
			system("ln --relative --symbolic $baiFile $bamFile.bai") if($samModule);
		} else {
			system("samtools index $bamFile");
		}
	}
}
{
	open(my $reader, ($variantFile =~ /\.gz$/ ? "gzip -dc $variantFile |" : $variantFile));
	chomp(my $line = <$reader>);
	$line =~ s/^#//;
	my @columnList = split(/\t/, $line, -1);
	print '#', join("\t", @columnList, 'Depth', 'Reference base depth', 'Alternate base depth', 'Variant allele frequency'), "\n";
	my @tokenListList = ();
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		push(@tokenListList, \@tokenList);
		if(scalar(@tokenListList) >= $numberPerThread) {
			if($threads == 1) {
				printTokenListList(@tokenListList);
			} else {
				forkPrintSubroutine(\&printTokenListList, @tokenListList);
			}
			@tokenListList = ();
		}
	}
	if(@tokenListList) {
		if($threads == 1) {
			printTokenListList(@tokenListList);
		} else {
			forkPrintSubroutine(\&printTokenListList, @tokenListList);
		}
		@tokenListList = ();
	}
	forkPrintWait();
	close($reader);

	sub printTokenListList {
		my $sam = Bio::DB::Sam->new(-bam => $bamFile) if($samModule);
		foreach(@_) {
			my @tokenList = @$_;
			my %tokenHash = ();
			@tokenHash{@columnList} = @tokenList;
			my ($chromosome, $position, $refBase, $altBase) = @tokenHash{'Chromosome', 'Position', 'Reference base', 'Alternate base'};
			($altBase, my $allele) = getAltBaseAllele($altBase);
			if($altBase eq '*') {
				$altBase = '';
			} else {
				($chromosome, $position, $refBase, $altBase) = extendIndel($chromosome, $position, $refBase, $altBase);
			}
			my ($start, $end, $strand) = ($position, $position + length($refBase) - 1, $tokenHash{'Strand'});
			my @readBaseList = ();
			if($samModule) {
				foreach my $alignment ($sam->get_features_by_location(-seq_id => $chromosome, -start => $end, -end => $end)) {
					next if($alignment->qual < $minimumMappingQuality);
					next if(($alignment->flag & int($includeFlag)) != $includeFlag);
					next if($alignment->flag & int($excludeFlag));
					next if($alignment->start > $start);
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
						next unless(grep {($alignment->flag & 253) == $_} 97, 145, 0);
					}
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
						next unless(grep {($alignment->flag & 253) == $_} 81, 161, 16);
					}
					my @positionList = getPositionList($alignment->start, $alignment->cigar_str);
					my $index = max(0, map {$_ + 1} grep {$positionList[$_] ne '' && $positionList[$_] < $start} 0 .. $#positionList);
					my $length = min(scalar(@positionList), grep {$positionList[$_] ne '' && $positionList[$_] > $end} 0 .. $#positionList) - $index;
					push(@readBaseList, substr($alignment->query->dna, $index, $length)) if($length >= 0);
				}
			} else {
				open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$end-$end |");
				while(my $line = <$reader>) {
					chomp($line);
					my %tokenHash = ();
					(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
					$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
					next if($tokenHash{'pos'} > $start);
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
						next unless(grep {($tokenHash{'flag'} & 253) == $_} 97, 145, 0);
					}
					if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
						next unless(grep {($tokenHash{'flag'} & 253) == $_} 81, 161, 16);
					}
					my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
					my $index = max(0, map {$_ + 1} grep {$positionList[$_] ne '' && $positionList[$_] < $start} 0 .. $#positionList);
					my $length = min(scalar(@positionList), grep {$positionList[$_] ne '' && $positionList[$_] > $end} 0 .. $#positionList) - $index;
					push(@readBaseList, substr($tokenHash{'seq'}, $index, $length)) if($length >= 0);
				}
				close($reader);
			}
			my $depth = scalar(@readBaseList);
			my $refBaseDepth = scalar(grep {$_ eq $refBase} @readBaseList);
			my $altBaseDepth = scalar(grep {$_ eq $altBase} @readBaseList);
			forkPrint(join("\t", @tokenList, $depth, $refBaseDepth, $altBaseDepth, $depth > 0 ? $altBaseDepth / $depth : 'NaN'), "\n");
		}
	}
}
if($copyBam) {
	system("rm -f $temporaryPrefix.bam");
	system("rm -f $temporaryPrefix.bam.bai");
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		}
	}
	return @positionList;
}

sub extendIndel {
	my ($chromosome, $position, $refBase, $altBase) = @_;
	if($refBase ne $altBase) {
		my $isDeletion = $refBase =~ /^$altBase/ || $refBase =~ /$altBase$/;
		while($refBase =~ /^$altBase/ || $altBase =~ /^$refBase/) {
			last if($position + length($refBase) - 1 == $db->length($chromosome));
			my $extBase = uc($db->seq($chromosome, ($_ = $position + length($refBase)), $_));
			($refBase, $altBase) = map {"$_$extBase"} ($refBase, $altBase);
		}
		while($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/) {
			last if($position == 1);
			my $extBase = uc($db->seq($chromosome, ($position = $position - 1), $position));
			($refBase, $altBase) = map {"$extBase$_"} ($refBase, $altBase);
		}
		if($isDeletion) {
			$position += 1;
			substr($_, -1, 1, '') foreach($refBase, $altBase);
			substr($_,  0, 1, '') foreach($refBase, $altBase);
		}
	}
	return ($chromosome, $position, $refBase, $altBase);
}

sub getAltBaseAllele {
	my ($altBase) = @_;
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	return ($altBase, 1) if(scalar(@altBaseList) == 1);
	foreach my $altBaseIndex (0 .. $#altBaseList) {
		return ($1, $altBaseIndex + 1) if($altBaseList[$altBaseIndex] =~ /^\[(.*)\]$/);
	}
}
