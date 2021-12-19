# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max);
use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $header = ''),
	'n' => \(my $noAlleleColumn = ''),
	'p' => \(my $positionOnly = ''),
	'a' => \(my $annotatedOnly = ''),
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
	my @columnList = getTokenList();
	my %columnIndexHash = getColumnIndexHash(@columnList);
	my @tokenList = getTokenList();
	my @tokenListList = ();
	sub getTokenListList {
		my ($chromosome, $start, $end) = @_;
		@tokenListList = grep {$_->[$columnIndexHash{'chromosome'}] eq $chromosome && $start <= $_->[$columnIndexHash{'end'}]} @tokenListList;
		while(@tokenList && defined($chromosomeIndexHash{$chromosome}) && $chromosomeIndexHash{$tokenList[$columnIndexHash{'chromosome'}]} < $chromosomeIndexHash{$chromosome}) {
			@tokenList = getTokenList();
		}
		while(@tokenList && $tokenList[$columnIndexHash{'chromosome'}] eq $chromosome && $tokenList[$columnIndexHash{'start'}] <= $end) {
			push(@tokenListList, [@tokenList]) if($start <= $tokenList[$columnIndexHash{'end'}]);
			@tokenList = getTokenList();
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
my %columnIndexHash = ();
my $vcf = '';
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	if($header ne '') {
		print join("\t", $line, @columnList), "\n";
		{
			$line =~ s/^#//;
			my @columnList = split(/\t/, $line, -1);
			%columnIndexHash = getColumnIndexHash(@columnList);
		}
		@columnList = grep {$_ ne 'allele'} @columnList if(defined($columnIndexHash{'allele'}));
		$header = '';
		next;
	}
	if($line =~ /^#/) {
		print "$line\n";
		if($line =~ /^##fileformat=VCF/) {
			$vcf = 1;
			$title = shift(@columnList) if(@columnList);
		}
		if($line =~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO/) {
			my @columnList = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO');
			%columnIndexHash = getColumnIndexHash(@columnList);
		}
		next;
	}
	my @tokenList = split(/\t/, $line, -1);
	my ($chromosome, $position, $refBase, $altBase) = @tokenList[%columnIndexHash ? @columnIndexHash{'chromosome', 'position', 'refBase', 'altBase'} : $vcf ? (0, 1, 3, 4) : (0, 1, 2, 3)];
	($refBase, $altBase) = map {uc} ($refBase, $altBase);
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	my @altBaseAlleleList = map {[$altBaseList[$_], $_ + 1]} 0 .. $#altBaseList;
	if(my @indexList = grep {$altBaseAlleleList[$_]->[0] =~ s/^\[(.*)\]$/$1/} 0 .. $#altBaseAlleleList) {
		@altBaseAlleleList = @altBaseAlleleList[@indexList];
	}
	my @tokenHashList = ();
	foreach(@altBaseAlleleList) {
		my ($altBase, $allele) = @$_;
		if($vcf eq '' && defined(my $index = $columnIndexHash{'allele'})) {
			next if($tokenList[$index] ne $allele);
		}
		if($altBase =~ /^[A-Za-z]*$/) {
			my ($start, $end) = getExtendedStartEnd($chromosome, stripIdentical(0, $position, $refBase, $altBase));
			foreach my $tokenList (getTokenListList($chromosome, $position, $end)) {
				my %tokenHash = ();
				$tokenHash{'allele'} = $allele if(scalar(@altBaseList) > 1);
				@tokenHash{(getColumnList())} = @$tokenList;
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
	my @columnList = ('allele', @columnList) if($noAlleleColumn eq '');
	if($vcf) {
		if(@tokenHashList) {
			if($title ne '') {
				my @infoList = grep {$_ ne '.' && $_ !~ /^$title\./} split(/;/, $tokenList[$columnIndexHash{'INFO'}]);
				for(my $number = 1; $number <= scalar(@tokenHashList); $number++) {
					my %tokenHash = %{$tokenHashList[$number - 1]};
					if(scalar(@tokenHashList) > 1) {
						push(@infoList, map {"$title.${_}_$number=$tokenHash{$_}"} grep {defined($tokenHash{$_}) && $tokenHash{$_} ne ''} @columnList);
					} else {
						push(@infoList, map {"$title.$_=$tokenHash{$_}"} grep {defined($tokenHash{$_}) && $tokenHash{$_} ne ''} @columnList);
					}
				}
				$tokenList[$columnIndexHash{'INFO'}] = join(';', @infoList);
			}
			print join("\t", @tokenList), "\n";
		} elsif($annotatedOnly eq '') {
			print join("\t", @tokenList), "\n";
		}
	} else {
		if(@tokenHashList) {
			open(my $writer, "| sort -u");
			print $writer join("\t", @tokenList, map {defined($_) ? $_ : ''} @$_{@columnList}), "\n" foreach(@tokenHashList);
			close($writer);
		} elsif($annotatedOnly eq '') {
			print join("\t", @tokenList, ('' x scalar(@columnList))), "\n";
		}
	}
}
close($reader);
closeAnnotationTableFileReader;

sub getColumnIndexHash {
	my (@columnList) = @_;
	my %columnIndexHash = ();
	@columnIndexHash{@columnList} = 0 .. $#columnList;
	foreach my $column ('chromosome', 'Chromosome', 'CHROM', 'chr') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'chromosome'} = $index;
			last;
		}
	}
	foreach my $column ('position', 'Position', 'POS') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'position'} = $index;
			last;
		}
	}
	foreach my $column ('refBase', 'haplotypeReference', 'Reference base', 'REF') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'refBase'} = $index;
			last;
		}
	}
	foreach my $column ('altBase', 'haplotypeAlternate', 'Alternate base', 'ALT') {
		if(defined(my $index = $columnIndexHash{$column})) {
			$columnIndexHash{'altBase'} = $index;
			last;
		}
	}
	return %columnIndexHash;
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
