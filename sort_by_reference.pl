# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $header = ''),
	'a' => \(my $append = ''),
);
if(scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl sort_by_reference.pl [options] input.txt reference.fasta chromosome.index position.index [...] > output.sorted.txt

Options: -h       the first line is header
         -a       append to files instead of open file writers

EOF
}
my ($tableFile, $referenceFastaFile, $chromosomeIndex, @positionIndexList) = @ARGV;
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
chomp(my $hostname = `hostname`);
system(sprintf('rm -rf %s', getChromosomeTemporaryFile('*')));
my %chromosomeWriterHash = ();
open(my $reader, $tableFile);
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^#/) {
		print "$line\n";
		next;
	}
	if($header) {
		$header = '';
		print "$line\n";
		next;
	}
	my @tokenList = split("\t", $line, -1);
	my $chromosome = $tokenList[$chromosomeIndex];
	$chromosome =~ s/[ |>;()\$].*$//;
	if(defined($chromosomeWriterHash{$chromosome})) {
		print {$chromosomeWriterHash{$chromosome}} "$line\n";
	} else {
		open(my $writer, ($append ? '>>' : '>'), getChromosomeTemporaryFile($chromosome));
		print $writer "$line\n";
		if($append) {
			close($writer);
		} else {
			$chromosomeWriterHash{$chromosome} = $writer;
		}
	}
}
close($reader);
foreach my $chromosome (keys %chromosomeWriterHash) {
	close($chromosomeWriterHash{$chromosome});
}
my @chromosomeList = getChromosomeList();
s/[ |>;()\$].*$// foreach(@chromosomeList);
my @sortOptionList = map {sprintf('-k%d,%dn', $_ + 1, $_ + 1)} @positionIndexList;
@sortOptionList = (sprintf('-k%d,%d', $chromosomeIndex + 1, $chromosomeIndex + 1), @sortOptionList);
foreach my $chromosome (@chromosomeList) {
	system("sort --field-separator='\t' @sortOptionList " . getChromosomeTemporaryFile($chromosome)) if(-e getChromosomeTemporaryFile($chromosome));
}
system(sprintf('rm -rf %s', getChromosomeTemporaryFile('*')));

sub getChromosomeTemporaryFile {
	my ($chromosome) = @_;
	return "$temporaryDirectory/sort_by_reference.$hostname.$$.$chromosome";
}

sub getChromosomeList {
	my @chromosomeList = ();
	if(my $faiFile = $referenceFastaFile =~ /\.fai$/ ? $referenceFastaFile : `find $referenceFastaFile.fai -newer $referenceFastaFile 2> /dev/null`) {
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
