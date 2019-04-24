# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long;

GetOptions(
	'c' => \(my $comment = ''),
	'h' => \(my $header = ''),
	'a' => \(my $append = ''),
);
my ($tableFile, $referenceFastaFile, $chromosomeIndex, @positionIndexList) = @ARGV;
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
chomp(my $hostname = `hostname`);
system(sprintf('rm -rf %s', getChromosomeTemporaryFile('*')));
my %chromosomeWriterHash = ();
open(my $reader, $tableFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^#/ && (!$comment || (print "$line\n")));
	next if($header && !($header = 0) && (print "$line\n"));
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
