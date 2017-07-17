# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;

my ($gtfFile, $geneInfoFile, $geneRefseqFile, $taxId) = @ARGV;
my %geneNameHash = ();
{
	open(my $reader, ($geneInfoFile =~ /\.gz$/) ? "gzip -dc $geneInfoFile |" : $geneInfoFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my @tokenList = split(/\t/, $line);
		if($tokenList[0] == $taxId) {
			$geneNameHash{$tokenList[1]} = $tokenList[2];
		}
	}
	close($reader);
}
my %transcriptGeneNameHash = ();
my %transcriptProteinHash = ();
{
	open(my $reader, ($geneRefseqFile =~ /\.gz$/) ? "gzip -dc $geneRefseqFile |" : $geneRefseqFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my @tokenList = split(/\t/, $line);
		if($tokenList[0] == $taxId) {
			(my $transcriptId = $tokenList[3]) =~ s/\.[0-9]+$//;
			(my $proteinId = $tokenList[5]) =~ s/\.[0-9]+$//;
			if($transcriptId ne '-') {
				$transcriptGeneNameHash{$transcriptId} = $geneNameHash{$tokenList[1]};
				$transcriptProteinHash{$transcriptId} = $proteinId if($proteinId ne '-');
			}
		}
	}
	close($reader);
}
open(my $reader, ($gtfFile =~ /\.gz$/) ? "gzip -dc $gtfFile |" : $gtfFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/);
	my %tokenHash = ();
	@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
	my @attributeList = ();
	my %attributeHash = ();
	while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g) {
		push(@attributeList, $1);
		$attributeHash{$1} = $2;
	}
	my $transcriptId = $attributeHash{'transcript_id'};
	$transcriptId =~ s/_dup[0-9]+$//;
	$transcriptId =~ s/\.[0-9]+$//;
	unless(defined($attributeHash{'gene_name'})) {
		if(defined($attributeHash{'gene_name'} = $transcriptGeneNameHash{$transcriptId})) {
			push(@attributeList, 'gene_name');
		}
	}
	unless(defined($attributeHash{'protein_id'})) {
		if(defined($attributeHash{'protein_id'} = $transcriptProteinHash{$transcriptId})) {
			push(@attributeList, 'protein_id');
		}
	}
	$tokenHash{'attribute'} = join(' ', map {"$_ \"$attributeHash{$_}\";"} @attributeList);
	print join("\t", @tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'}), "\n";
}
