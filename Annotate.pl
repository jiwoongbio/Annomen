# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max);
use Bio::DB::Fasta;
use Getopt::Long;

GetOptions('p' => \(my $positionOnly = ''));
my ($inputFile, $referenceFastaFile, $annotationTableFile, @columnList) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
{
	chomp(my @chromosomeList = `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null | xargs cat | cut -f1`);
	chomp(@chromosomeList = `grep '^>' $referenceFastaFile | sed 's/^>//' | sed 's/\\s.*\$//'`) unless(@chromosomeList);
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
	$columnIndexHash{'chromosome'} = $columnIndexHash{'chr'} unless(defined($columnIndexHash{'chromosome'}));
	my @tokenList = getTokenList;
	my @tokenListList = ();
	sub setChromosomeStartEnd {
		my ($chromosome, $start, $end) = @_;
		@tokenListList = grep {$_->[$columnIndexHash{'chromosome'}] eq $chromosome && $start <= $_->[$columnIndexHash{'end'}]} @tokenListList;
		while(@tokenList && defined($chromosomeIndexHash{$chromosome}) && $chromosomeIndexHash{$tokenList[$columnIndexHash{'chromosome'}]} < $chromosomeIndexHash{$chromosome}) {
			@tokenList = getTokenList;
		}
		while(@tokenList && $tokenList[$columnIndexHash{'chromosome'}] eq $chromosome && $tokenList[$columnIndexHash{'start'}] <= $end) {
			push(@tokenListList, [@tokenList]) if($start <= $tokenList[$columnIndexHash{'end'}]);
			@tokenList = getTokenList;
		}
	}
	sub getTokenListList {
		if($positionOnly) {
			my ($start, $end) = @_;
			return grep {$_->[$columnIndexHash{'start'}] <= $end} @tokenListList;
		} else {
			my ($position, $refBase, $altBase) = @_;
			return grep {$_->[$columnIndexHash{'start'}] == $position && $_->[$columnIndexHash{'haplotypeReference'}] eq $refBase && $_->[$columnIndexHash{'haplotypeAlternate'}] eq $altBase} @tokenListList;
		}
	}
	sub getColumnList {
		return @columnList;
	}
	sub closeAnnotationTableFileReader {
		close($reader);
	}
}
my ($vcf, $title) = (0);
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^#/) {
		$title = shift(@columnList) if($line =~ /^##fileformat=VCF/ && ($vcf = 1));
		print "$line\n";
	} else {
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $position, $refBase, $altBase) = @tokenList[$vcf ? (0, 1, 3, 4) : (0, 1, 2, 3)];
		my @tokenListList = ();
		if($positionOnly) {
			my @startEndList = map {[getExtendedStartEnd($chromosome, stripIdentical(0, $position, $refBase, $_))]} split(/,/, $altBase);
			setChromosomeStartEnd($chromosome, $position, my $end = max(map {$_->[1]} @startEndList));
			push(@tokenListList, getTokenListList(min(map {$_->[0]} @startEndList), $end));
		} else {
			setChromosomeStartEnd($chromosome, $position, $position + length($refBase) - 1);
			push(@tokenListList, getTokenListList(stripIdentical(1, $position, $refBase, $_))) foreach(split(/,/, $altBase));
		}
		my @tokensList = ();
		foreach my $tokenList (@tokenListList) {
			my %tokenHash = ();
			@tokenHash{(getColumnList)} = @$tokenList;
			push(@tokensList, join("\t", map {defined($_) ? $_ : ''} @tokenHash{@columnList}));
		}
		@tokensList = sort @tokensList;
		@tokensList = map {$tokensList[$_]} grep {$_ == 0 || $tokensList[$_ - 1] ne $tokensList[$_]} 0 .. $#tokensList;
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
	my ($head, $position, @sequenceList) = @_;
	while(scalar(grep {length($_) <= $head} @sequenceList) == 0 && (my @baseList = map {substr($_, -1, 1)} @sequenceList)) {
		last if(grep {$baseList[$_ - 1] ne $baseList[$_]} 1 .. $#baseList);
		substr($_, -1, 1, '') foreach(@sequenceList);
	}
	while(scalar(grep {length($_) <= $head} @sequenceList) == 0 && (my @baseList = map {substr($_, 0, 1)} @sequenceList)) {
		last if(grep {$baseList[$_ - 1] ne $baseList[$_]} 1 .. $#baseList);
		substr($_, 0, 1, '') foreach(@sequenceList);
		$position += 1;
	}
	return ($position, @sequenceList);
}
