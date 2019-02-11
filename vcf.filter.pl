# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use Getopt::Long;

GetOptions('g' => \(my $checkGenotypes = ''), 'v' => \(my $checkVariantsOnly = ''), 'r' => \(my $reverse = ''), 'd=s' => \(my $defaults = ''));
my ($vcfFile, $filterName, @filterExpressionList) = @ARGV;
for(my $index = 0; $index < scalar(@filterExpressionList); $index++) {
	$filterExpressionList[$index] =~ s/([A-Za-z_.0-9]*[A-Za-z]+[A-Za-z_.0-9]*)/\$keyValueHash{'$1'}/g;
}
open(my $reader, $vcfFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line !~ /^#/) {
		my ($chromosome, $position, $id, $refBase, $altBase, $quality, $filter, $info, $format, @genotypeList) = split(/\t/, $line);
		my %filterNameHash = map {$_ ne '.' && $_ ne 'PASS' ? ($_ => 1) : ()} split(/;/, $filter);
		$filterNameHash{$filterName} = 0;
		if($checkGenotypes) {
			my $count = 0;
			foreach my $genotype (@genotypeList) {
				my %keyValueHash = ();
				$keyValueHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} split(/;/, $defaults));
				@keyValueHash{split(/:/, $format)} = split(/:/, $genotype);
				$keyValueHash{'DP'} = 0 if(!defined($keyValueHash{'DP'}) || $keyValueHash{'DP'} eq '.');
				$keyValueHash{'GQ'} = 0 if(!defined($keyValueHash{'GQ'}) || $keyValueHash{'GQ'} eq '.');
				next if($checkVariantsOnly && $keyValueHash{'GT'} eq '0/0');
				$count++ if(map {eval($_) ? $_ : ()} @filterExpressionList);
			}
			if($reverse) {
				$filterNameHash{$filterName} = 1 if($count == 0);
			} else {
				$filterNameHash{$filterName} = 1 if($count > 0);
			}
		} elsif(not ($checkVariantsOnly && $altBase eq '.')) {
			my %keyValueHash = ('QUAL' => $quality);
			foreach my $keyValue (split(/;/, "$defaults;$info")) {
				$keyValueHash{$1} = $2 if($keyValue =~ /^(.+)=(.+)$/);
			}
			$keyValueHash{'QD'} = 0 if(!defined($keyValueHash{'QD'}) || $keyValueHash{'QD'} eq '.');
			$filterNameHash{$filterName} = 1 if(map {eval($_) ? $_ : ()} @filterExpressionList);
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
