# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use Getopt::Long qw(:config no_ignore_case);

my @sampleList = ();
GetOptions(
	'h' => \(my $help = ''),
	's=s' => \@sampleList,
	'H' => \(my $hasHeader = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl variant_type_count.pl [options] variant.txt [...] > variant_type_count.txt

Options: -h       display this help message
         -s STR   sample
         -H       has header

EOF
}
@sampleList = map {split(/,/, $_)} @sampleList;
my %sampleHash = map {$_ => 1} @sampleList;

my %typeExpressionHash = ();
$typeExpressionHash{'total'} = '1';

$typeExpressionHash{'homozygous'}   = 'scalar(getHaplotypeList($genotype)) == 1';
$typeExpressionHash{'heterozygous'} = 'scalar(getHaplotypeList($genotype)) > 1';

$typeExpressionHash{'transition'}   = 'grep {$_->[0] eq $refBase && $_->[1] eq $altBase} ["A", "G"], ["C", "T"], ["G", "A"], ["T", "C"]';
$typeExpressionHash{'transversion'} = 'grep {$_->[0] eq $refBase && $_->[1] eq $altBase} ["A", "C"], ["A", "T"], ["C", "A"], ["C", "G"], ["G", "C"], ["G", "T"], ["T", "A"], ["T", "G"]';

my @regionList = ('CDS', 'utr5', 'utr3', 'intron', 'non_coding_exon', 'non_coding_intron');
$typeExpressionHash{$_} = "\$region =~ /$_/" foreach(@regionList);
$typeExpressionHash{'exon'} = join(' || ', map {$typeExpressionHash{$_}} ('CDS', 'utr5', 'utr3', 'non_coding_exon'));

my @mutationList = ('silent', 'missense', 'inframe', 'frameshift', 'nonsense', 'readthrough', 'startcodon', 'splicing', 'junction');
$typeExpressionHash{$_} = "\$mutation =~ /$_/" foreach(@mutationList);
$typeExpressionHash{'NS/SS/I'} = join(' || ', map {$typeExpressionHash{$_}} grep {$_ ne 'silent'} @mutationList);

my @typeList = ('total', 'homozygous', 'heterozygous', 'transition', 'transversion', 'exon', @regionList, @mutationList, 'NS/SS/I');

my (@variantFileList) = @ARGV;
my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1 -k2,2 -k3,6 | uniq | cut -f1,2 | uniq -c");
foreach my $variantFile (@variantFileList) {
	my %columnIndexHash = ();
	open(my $reader, $variantFile);
	if($hasHeader) {
		chomp(my $line = <$reader>);
		$line =~ s/^#//;
		my @columnList = split(/\t/, $line, -1);
		@columnIndexHash{@columnList} = 0 .. $#columnList;
		($columnIndexHash{'chromosome'}) = grep {defined} @columnIndexHash{'chromosome', 'Chromosome'};
		($columnIndexHash{'position'}) = grep {defined} @columnIndexHash{'position', 'Position'};
		($columnIndexHash{'refBase'}) = grep {defined} @columnIndexHash{'refBase', 'Reference base'};
		($columnIndexHash{'altBase'}) = grep {defined} @columnIndexHash{'altBase', 'Alternate base'};
		($columnIndexHash{'sample'}) = grep {defined} @columnIndexHash{'sample', 'Sample'};
		($columnIndexHash{'genotype'}) = grep {defined} @columnIndexHash{'genotype', 'Genotype', 'Tumor genotype'};
		($columnIndexHash{'gene'}) = grep {defined} @columnIndexHash{'gene', 'Gene name'};
		($columnIndexHash{'region'}) = grep {defined} @columnIndexHash{'region', 'Region type'};
		($columnIndexHash{'mutation'}) = grep {defined} @columnIndexHash{'mutation', 'Mutation class'};
	} else {
		my @columnList = ('chromosome', 'position', 'refBase', 'altBase', 'sample', 'genotype', 'gene', 'region', 'mutation');
		@columnIndexHash{@columnList} = 0 .. $#columnList;
	}
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my @tokenList = split(/\t/, $line, -1);
		my ($chromosome, $position, $refBase, $altBase, $sample, $genotype, $gene, $region, $mutation) = map {defined($_) ? $tokenList[$_] : $_} @columnIndexHash{'chromosome', 'position', 'refBase', 'altBase', 'sample', 'genotype', 'gene', 'region', 'mutation'};
		if(defined($sample)) {
			if(@sampleList) {
				next unless($sampleHash{$sample});
			}
		} else {
			if(scalar(@sampleList) == 1) {
				$sample = $sampleList[0];
			}
		}
		foreach my $type (@typeList) {
			my $expression = $typeExpressionHash{$type};
			print $writer join("\t", $sample, $type, $chromosome, $position, $refBase, $altBase), "\n" if(eval($expression));
		}
	}
	close($reader);
}
close($writer);
my %sampleTypeCountHash = ();
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ /^ *([0-9]+) (.*)$/) {
		my $count = $1;
		my ($sample, $type) = split(/\t/, $2, -1);
		$sampleTypeCountHash{$sample}->{$type} = $count;
	}
}
close($reader);
waitpid($pid, 0);
@sampleList = sort keys %sampleTypeCountHash unless(@sampleList);
print join("\t", 'sample', @sampleList), "\n";
foreach my $type (@typeList) {
	print join("\t", $type, map {defined($_) ? $_ : 0} map {$_->{$type}} @sampleTypeCountHash{@sampleList}), "\n";
}

sub getHaplotypeList {
	my ($genotype) = @_;
	my %haplotypeHash = ();
	$haplotypeHash{$_} = 1 foreach(split(/[\/|]/, $genotype));
	my @haplotypeList = sort {$a <=> $b} grep {$_ ne '.'} keys %haplotypeHash;
	return @haplotypeList;
}
