# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min);
use Bio::DB::Fasta;

my ($inputFile, $referenceFastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);
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
		($refBase, $altBase) = map {uc} ($refBase, $altBase);
		my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase);
		($position, $refBase, @altBaseList) = leftalignIndel($chromosome, stripIdentical($position, $refBase, @altBaseList));
		@tokenList[$vcf ? (0, 1, 3, 4) : (0, 1, 2, 3)] = ($chromosome, $position, $refBase, $altBase = join(',', @altBaseList));
		print join("\t", @tokenList), "\n";
	}
}
close($reader);

sub leftalignIndel {
	my ($chromosome, $position, @sequenceList) = @_;
	if(grep {$_ eq ''} @sequenceList) {
		while(my $base = uc($db->seq($chromosome, $position = $position - 1, $position))) {
			@sequenceList = map {"$base$_"} @sequenceList;
			my @baseList = map {substr($_, -1, 1)} @sequenceList;
			last if(grep {$baseList[$_ - 1] ne $baseList[$_]} 1 .. $#baseList);
			substr($_, -1, 1, '') foreach(@sequenceList);
		}
	}
	return ($position, @sequenceList);
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
