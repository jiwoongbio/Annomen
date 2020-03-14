# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max);
use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'p' => \(my $positionOnly = ''),
);
my ($inputFile, $referenceFastaFile, $annotationTableFile, @columnList) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
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
	$columnIndexHash{'chromosome'} = $columnIndexHash{'chr'} unless(defined($columnIndexHash{'chromosome'}));
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
my $title = '';
@columnList = ('allele', @columnList);
my $vcf = '';
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^#/) {
		if($line =~ /^##fileformat=VCF/) {
			$vcf = 1;
			$title = splice(@columnList, 1, 1);
		}
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
			my ($start, $end) = getExtendedStartEnd($chromosome, stripIdentical(0, $position, $refBase, $altBase));
			foreach my $tokenList (getTokenListList($chromosome, $position, $end)) {
				my %tokenHash = ();
				$tokenHash{'allele'} = $altBaseIndex + 1 if(scalar(@altBaseList) > 1);
				@tokenHash{(getColumnList)} = @$tokenList;
				next if($tokenHash{'end'} < $start);
				if($positionOnly) {
					push(@tokenHashList, \%tokenHash);
				} else {
					my ($position, $refBase, $altBase) = stripIdentical(1, $position, $refBase, $altBase);
					push(@tokenHashList, \%tokenHash) if($tokenHash{'start'} == $position && $tokenHash{'haplotypeReference'} eq $refBase && $tokenHash{'haplotypeAlternate'} eq $altBase);
				}
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
