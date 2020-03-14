# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my @defaultList = ();
GetOptions(
	'h' => \(my $help = ''),
	'g' => \(my $filterGenotype = ''),
	'v' => \(my $invert = ''),
	'd=s' => \@defaultList,
);
my %defaultTokenHash = ();
$defaultTokenHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_, 2)]} @defaultList);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vcf.filter.pl [options] input.vcf filter.name filter.expression [...]

Options: -h       display this help message
         -g       filter genotype field
         -v       invert
         -d STR   degault value, key=value

EOF
}
my ($vcfFile, $filterName, @filterExpressionList) = @ARGV;
my %filterTokenHash = ();
foreach(@filterExpressionList) {
	$filterTokenHash{$1} = 1 while(/([A-Za-z_.0-9]*[A-Za-z]+[A-Za-z_.0-9]*)/g);
	s/([A-Za-z_.0-9]*[A-Za-z]+[A-Za-z_.0-9]*)/\$tokenHash{'$1'}/g;
}
open(my $reader, ($vcfFile =~ /\.gz$/ ? "gzip -dc $vcfFile |" : $vcfFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line !~ /^#/) {
		my ($chromosome, $position, $id, $refBase, $altBase, $quality, $filter, $info, $format, @genotypeList) = split(/\t/, $line, -1);
		my %filterNameHash = map {$_ => 1} grep {$_ ne '.' && $_ ne 'PASS'} split(/;/, $filter, -1);
		$filterNameHash{$filterName} = '';
		my $filterCount = 0;
		if($filterGenotype eq '') {
			my %tokenHash = %defaultTokenHash;
			$tokenHash{'QUAL'} = $quality;
			$tokenHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_, 2)]} split(/;/, $info, -1));
			$tokenHash{$_} = 0 foreach(grep {!defined($tokenHash{$_}) || $tokenHash{$_} eq '.'} keys %filterTokenHash);
			$filterCount += 1 if(grep {eval($_)} @filterExpressionList);
		} else {
			foreach my $genotype (@genotypeList) {
				my %tokenHash = %defaultTokenHash;
				@tokenHash{split(/:/, $format, -1)} = split(/:/, $genotype, -1);
				$tokenHash{$_} = 0 foreach(grep {!defined($tokenHash{$_}) || $tokenHash{$_} eq '.'} keys %filterTokenHash);
				$filterCount += 1 if(grep {eval($_)} @filterExpressionList);
			}
		}
		if($invert eq '') {
			$filterNameHash{$filterName} = 1 if($filterCount > 0);
		} else {
			$filterNameHash{$filterName} = 1 if($filterCount == 0);
		}
		my @filterNameList = grep {$filterNameHash{$_}} sort keys %filterNameHash;
		$filter = @filterNameList ? join(';', @filterNameList) : 'PASS';
		if(defined($format)) {
			$line = join("\t", $chromosome, $position, $id, $refBase, $altBase, $quality, $filter, $info, $format, @genotypeList);
		} else {
			$line = join("\t", $chromosome, $position, $id, $refBase, $altBase, $quality, $filter, $info);
		}
	}
	print "$line\n";
}
close($reader);
