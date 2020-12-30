# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

my @codonList = ();
GetOptions(
	'h' => \(my $help = ''),
	'C=s' => \@codonList,
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl NS_siteCount.pl [options] annotation_table.txt transcript.fasta > NS_siteCount.txt

Options: -h       display this help message
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 1 (standard)]

EOF
}
my %codonHash = (
	'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
	'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
	'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
	'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

	'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
	'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
	'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
	'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

	'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
	'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
	'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
	'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

	'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
	'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
	'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
	'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
);
$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);

my ($annotationTableFile, $transcriptFastaFile) = @ARGV;

my %transcriptSequenceHash = ();
{
	my $transcriptId = '';
	open(my $reader, ($transcriptFastaFile =~ /\.gz$/ ? "gzip -dc $transcriptFastaFile |" : $transcriptFastaFile));
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^>(\S*)/ && ($transcriptId = $1));
		$transcriptSequenceHash{$transcriptId} .= $line;
	}
	close($reader);
}

my $pid = open2(my $reader, my $writer, 'sort -u');
{
	open(my $reader, $annotationTableFile);
	chomp(my $line = <$reader>);
	my @columnList = split(/\t/, $line, -1);
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);
		print $writer join("\t", @tokenHash{'geneName', 'transcriptId', 'proteinId', 'codingRegions', 'frame'}), "\n" if($tokenHash{'codingRegions'} ne '');
	}
	close($reader);
}
close($writer);
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{'geneName', 'transcriptId', 'proteinId', 'codingRegions', 'frame'} = split(/\t/, $line, -1);
	my @baseList = split(//, $transcriptSequenceHash{$tokenHash{'transcriptId'}});
	my $codingSequence = join('', @baseList[map {$_ - 1} eval($tokenHash{'codingRegions'})]);
	if(my $frame = $tokenHash{'frame'}) {
		$codingSequence = substr($codingSequence, $frame) if($frame > 0);
		$codingSequence = ('N' x -$frame) . $codingSequence if($frame < 0);
	}
	$tokenHash{'synonymousSiteCount'} = 0;
	$tokenHash{'nonsynonymousSiteCount'} = 0;
	for(my $index = 0; $index < length($codingSequence); $index += 3) {
		my $codon = substr($codingSequence, $index, 3);
		if(defined(my $aa = $codonHash{$codon})) {
			foreach my $index (0 .. 2) {
				my $base = substr($codon, $index, 1);
				foreach my $base (grep {$_ ne $base} 'A', 'C', 'G', 'T') {
					substr($codon, $index, 1, $base);
					if($codonHash{$codon} ne $aa) {
						$tokenHash{'nonsynonymousSiteCount'} += 1 / 3;
					} else {
						$tokenHash{'synonymousSiteCount'} += 1 / 3;
					}
				}
				substr($codon, $index, 1, $base);
			}
		}
	}
	print join("\t", @tokenHash{'geneName', 'transcriptId', 'proteinId', 'nonsynonymousSiteCount', 'synonymousSiteCount'}), "\n";
}
close($reader);
waitpid($pid, 0);
