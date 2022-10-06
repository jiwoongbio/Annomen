# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Bio::DB::Fasta;
use Getopt::Long;

GetOptions(
	'c' => \(my $printCount = ''),
);
my ($referenceFastaFile, @variantFileList) = @ARGV;
my $db = Bio::DB::Fasta->new($referenceFastaFile);

my %trinucleotidesCountHash = ();
my $totalCount = 0;
foreach my $variantFile (@variantFileList) {
	my $pid = open2(my $reader, my $writer, "sort -u");
	{
		open(my $reader, ($variantFile =~ /\.gz$/ ? "gzip -dc $variantFile |" : $variantFile));
		chomp(my $line = <$reader>);
		$line =~ s/^#//;
		my @columnList = split(/\t/, $line, -1);
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			my ($chromosome, $position, $refBase, $altBase) = @tokenHash{'Chromosome', 'Position', 'Reference base', 'Alternate base'};
			($altBase, my $allele) = getAltBaseAllele($altBase);
			next unless(grep {$_ eq $allele} split(/[\/|]/, $tokenHash{'Genotype'}));
			if($refBase =~ /^[ACGT]$/ && $altBase =~ /^[ACGT]$/ && $refBase ne $altBase) {
				print $writer join("\t", $chromosome, $position, $refBase, $altBase), "\n";
			}
		}
		close($reader);
	}
	close($writer);
	{
		while(my $line = <$reader>) {
			chomp($line);
			my ($chromosome, $position, $refBase, $altBase) = split(/\t/, $line, -1);
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
	waitpid($pid, 0);
}

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

sub getAltBaseAllele {
	my ($altBase) = @_;
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	return ($altBase, 1) if(scalar(@altBaseList) == 1);
	foreach my $altBaseIndex (0 .. $#altBaseList) {
		return ($1, $altBaseIndex + 1) if($altBaseList[$altBaseIndex] =~ /^\[(.*)\]$/);
	}
}
