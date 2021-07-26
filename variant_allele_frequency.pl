# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max);
use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'q=i' => \(my $minimumMappingQuality = 0),
	'f=i' => \(my $includeFlag = 0),
	'F=i' => \(my $excludeFlag = 0),
	'S=s' => \(my $stranded = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl variant_allele_frequency.pl [options] variant.txt reference.fasta sample.bam > sample.variant_allele_frequency.txt

Options: -h       display this help message
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -f INT   include flag [$includeFlag]
         -F INT   exclude flag [$excludeFlag]
         -S STR   stranded, "forward", "f", "reverse", or "r"

EOF
}
my ($variantFile, $referenceFastaFile, $bamFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
{
	open(my $reader, ($variantFile =~ /\.gz$/ ? "gzip -dc $variantFile |" : $variantFile));
	chomp(my $line = <$reader>);
	$line =~ s/^#//;
	my @columnList = split(/\t/, $line, -1);
	print '#', join("\t", @columnList, 'Depth', 'Reference base depth', 'Alternate base depth', 'Variant allele frequency'), "\n";
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my %tokenHash = ();
		@tokenHash{@columnList} = @tokenList;
		my ($chromosome, $position, $refBase, $altBase) = @tokenHash{'Chromosome', 'Position', 'Reference base', 'Alternate base'};
		($altBase, my $allele) = getAltBaseAllele($altBase);
		if($altBase eq '*') {
			$altBase = '';
		} else {
			($chromosome, $position, $refBase, $altBase) = extendIndel($chromosome, $position, $refBase, $altBase);
		}
		my $strand = $tokenHash{'Strand'};
		my $depth = scalar(my @readBaseList = getReadBaseList($bamFile, $chromosome, $position, $position + length($refBase) - 1, $strand));
		my $refBaseDepth = scalar(grep {$_ eq $refBase} @readBaseList);
		my $altBaseDepth = scalar(grep {$_ eq $altBase} @readBaseList);
		print join("\t", @tokenList, $depth, $refBaseDepth, $altBaseDepth, $depth > 0 ? $altBaseDepth / $depth : 'NaN'), "\n";
	}
	close($reader);
}

sub getReadBaseList {
	my ($bamFile, $chromosome, $start, $end, $strand) = @_;
	my @readBaseList = ();
	open(my $reader, "samtools view -q $minimumMappingQuality -f $includeFlag -F $excludeFlag $bamFile $chromosome:$end-$end |");
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		next if($tokenHash{'pos'} > $start);
		if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '+') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '-')) {
			next unless(($tokenHash{'flag'} & 253) == 97 || ($tokenHash{'flag'} & 253) == 145 || ($tokenHash{'flag'} & 253) == 0);
		}
		if((($stranded eq 'f' || $stranded eq 'forward') && $strand eq '-') || (($stranded eq 'r' || $stranded eq 'reverse') && $strand eq '+')) {
			next unless(($tokenHash{'flag'} & 253) == 81 || ($tokenHash{'flag'} & 253) == 161 || ($tokenHash{'flag'} & 253) == 16);
		}
		my @positionList = getPositionList(@tokenHash{'pos', 'cigar'});
		my $index = max(0, map {$_ + 1} grep {$positionList[$_] ne '' && $positionList[$_] < $start} 0 .. $#positionList);
		my $length = min(scalar(@positionList), grep {$positionList[$_] ne '' && $positionList[$_] > $end} 0 .. $#positionList) - $index;
		push(@readBaseList, substr($tokenHash{'seq'}, $index, $length)) if($length >= 0);
	}
	close($reader);
	return @readBaseList;
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
			my $extBase = uc($db->seq($chromosome, ($_ = $position + length($refBase)), $_));
			($refBase, $altBase) = map {"$_$extBase"} ($refBase, $altBase);
		}
		while($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/) {
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
