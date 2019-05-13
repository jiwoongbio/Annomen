# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Bio::DB::Fasta;
use Getopt::Long;

GetOptions(
	'c' => \(my $printCount = ''),
);
my ($inputFile, $referenceFastaFile) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);

my %trinucleotidesCountHash = ();
my $totalCount = 0;
open(my $reader, $inputFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/);
	my ($chromosome, $position, $refBase, $altBase) = split(/\t/, $line, -1);
	($refBase, $altBase) = map {uc} ($refBase, $altBase);
	if($refBase =~ /^[ACGT]$/ && $altBase =~ /^[ACGT]$/ && $refBase ne $altBase) {
		my $refTrinucleotide = uc($db->seq($chromosome, $position - 1, $position + 1));
		substr(my $altTrinucleotide = $refTrinucleotide, 1, 1, $altBase);
		unless(grep {$_ eq $refBase} ('C', 'T')) {
			($refTrinucleotide = reverse($refTrinucleotide)) =~ tr/ACGT/TGCA/;
			($altTrinucleotide = reverse($altTrinucleotide)) =~ tr/ACGT/TGCA/;
		}
		my $trinucleotides = "$refTrinucleotide > $altTrinucleotide";
		$trinucleotidesCountHash{$trinucleotides} += 1;
		$totalCount += 1;
	}
}
close($reader);

foreach my $refBase ('C', 'T') {
	foreach my $altBase (grep {$_ ne $refBase} ('A', 'C', 'G', 'T')) {
		foreach my $base5prime ('A', 'C', 'G', 'T') {
			foreach my $base3prime ('A', 'C', 'G', 'T') {
				my $refTrinucleotide = "$base5prime$refBase$base3prime";
				my $altTrinucleotide = "$base5prime$altBase$base3prime";
				my $trinucleotides = "$refTrinucleotide > $altTrinucleotide";
				my $count = defined($_ = $trinucleotidesCountHash{$trinucleotides}) ? $_ : 0;
				print join("\t", "$base5prime\[$refBase>$altBase\]$base3prime", $printCount ? $count : $count / $totalCount), "\n";
			}
		}
	}
}
