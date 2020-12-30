# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl NS_mutationRatio.pl [options] NS_siteCount.txt variant.annotation.txt [...] > NS_mutationRatio.txt

Options: -h       display this help message

EOF
}
my ($siteCountFile, @variantFileList) = @ARGV;
my %mutationCountHash = ();
foreach my $variantFile (@variantFileList) {
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
		if(defined($tokenHash{'Genotype'})) {
			next unless(grep {$_ eq $allele} split(/[\/|]/, $tokenHash{'Genotype'}));
		}
		if(length($refBase) == 1 && length($altBase) == 1 && $tokenHash{'Region type'} eq 'CDS') {
			if($tokenHash{'Mutation class'} ne 'silent') {
				$mutationCountHash{join("\t", @tokenHash{'Gene name', 'Transcript ID', 'Protein ID'})}->{'nonsynonymous'} += 1;
			} else {
				$mutationCountHash{join("\t", @tokenHash{'Gene name', 'Transcript ID', 'Protein ID'})}->{'synonymous'} += 1;
			}
		}
	}
	close($reader);
}
{
	open(my $reader, $siteCountFile);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{'geneName', 'transcriptId', 'proteinId', 'nonsynonymousSiteCount', 'synonymousSiteCount'} = split(/\t/, $line, -1);
		if(defined(my $mutationCountHash = $mutationCountHash{join("\t", @tokenHash{'geneName', 'transcriptId', 'proteinId'})})) {
			@tokenHash{'nonsynonymousMutationCount', 'synonymousMutationCount'} = map {defined($_) ? $_ : 0} @$mutationCountHash{'nonsynonymous', 'synonymous'};
			my $dN_number = 1 - 4 / 3 * ($tokenHash{'nonsynonymousMutationCount'} / $tokenHash{'nonsynonymousSiteCount'});
			my $dN = $dN_number > 0 ? - 3 / 4 * log($dN_number) : 'NaN';
			my $dS_number = 1 - 4 / 3 * ($tokenHash{'synonymousMutationCount'} / $tokenHash{'synonymousSiteCount'});
			my $dS = $dS_number > 0 ? - 3 / 4 * log($dS_number) : 'NaN';
			print join("\t", @tokenHash{'geneName', 'transcriptId', 'proteinId', 'nonsynonymousSiteCount', 'synonymousSiteCount', 'nonsynonymousMutationCount', 'synonymousMutationCount'}, ($dN eq 'NaN' || $dS eq 'NaN') ? 'NaN' : $dS > 0 ? ($dN / $dS) : ($dN > 0 ? 'Inf' : 'NaN')), "\n";
		}
	}
	close($reader);
}

sub getAltBaseAllele {
	my ($altBase) = @_;
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	return ($altBase, 1) if(scalar(@altBaseList) == 1);
	foreach my $altBaseIndex (0 .. $#altBaseList) {
		return ($1, $altBaseIndex + 1) if($altBaseList[$altBaseIndex] =~ /^\[(.*)\]$/);
	}
}
