# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'v' => \(my $printVariantOnly = ''),
);
my ($vcfFile, @sampleList) = @ARGV;
my @columnList = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', @sampleList);
my %columnIndexHash = ();
open(my $reader, ($vcfFile =~ /\.gz$/ ? "gzip -dc $vcfFile |" : $vcfFile));
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^##/) {
		print "$line\n";
	} elsif($line =~ s/^#//) {
		my @tokenList = split(/\t/, $line);
		%columnIndexHash = map {$tokenList[$_] => $_} 0 .. $#tokenList;
		print '#', join("\t", @tokenList[@columnIndexHash{@columnList}]), "\n";
	} else {
		my @tokenList = split(/\t/, $line);
		if($printVariantOnly) {
			my @formatList = split(/:/, $tokenList[$columnIndexHash{'FORMAT'}], -1);
			my %formatIndexHash = map {$formatList[$_] => $_} 0 .. $#formatList;
			my @genotypeList = map {$_->[$formatIndexHash{'GT'}]} map {[split(/:/, $_)]} @tokenList[@columnIndexHash{@sampleList}];
			@genotypeList = unphaseGenotype(@genotypeList);
			print join("\t", @tokenList[@columnIndexHash{@columnList}]), "\n" if(grep {$_ ne '' && $_ ne '0/0'} @genotypeList);
		} else {
			print join("\t", @tokenList[@columnIndexHash{@columnList}]), "\n";
		}
	}
}
close($reader);

sub unphaseGenotype {
	return join('/', sort {$a <=> $b} grep {$_ ne '.'} split(/[\/|]/, $_[0])) if(scalar(@_) == 1);
	return map {unphaseGenotype($_)} @_;
}
