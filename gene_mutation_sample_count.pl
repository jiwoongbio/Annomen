# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my @mutationList = ('silent', 'missense', 'inframe', 'frameshift', 'nonsense', 'readthrough', 'startcodon', 'splicing', 'junction');

my @columnList = ();
my @additionalColumnValueFileList = ();
GetOptions(
	'h' => \(my $help = ''),
	'c=s' => \@columnList,
	'a=s' => \@additionalColumnValueFileList,
	's=s' => \(my $samplePatientFile = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl gene_mutation_sample_count.pl [options] group=variant.txt [...] > gene_mutation_sample_count.txt

Options: -h       display this help message
         -c STR   mutation column e.g. missense, NS/SS/I
         -a FILE  additional column-value file
         -s FILE  sample-patient file, counting patients instead of samples

EOF
}
push(@columnList, @mutationList, 'NS/SS/I') if(scalar(@columnList) == 0);
@additionalColumnValueFileList = map {[$_->[0], $_->[-1]]} map {[split(/=/, $_, 2)]} @additionalColumnValueFileList;
my %geneGroupColumnValueHash = ();
foreach(@additionalColumnValueFileList) {
	my ($column, $valueFile) = @$_;
	push(@columnList, $column);
	open(my $reader, ($valueFile =~ /\.gz$/ ? "gzip -dc $valueFile |" : $valueFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($group, $gene, $value) = split(/\t/, $line, -1);
		$geneGroupColumnValueHash{$gene}->{$group}->{$column} = $value;
	}
	close($reader);
}
my %samplePatientHash = ();
if($samplePatientFile ne '') {
	open(my $reader, ($samplePatientFile =~ /\.gz$/ ? "gzip -dc $samplePatientFile |" : $samplePatientFile));
	while(my $line = <$reader>) {
		chomp($line);
		my ($sample, $patient) = split(/\t/, $line, -1);
		$samplePatientHash{$sample} = $patient;
	}
	close($reader);
}
my @groupVariantFileList = ();
foreach(@ARGV) {
	my ($group, $variantFile) = split(/=/, $_, 2);
	chomp(my @variantFileList = `ls $variantFile`);
	push(@groupVariantFileList, [$group, $_]) foreach(@variantFileList);
}
my %groupHash = ();
foreach(@groupVariantFileList) {
	my ($group, $variantFile) = @$_;
	$groupHash{$group} = 1;
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
		my $variant = join("\t", $chromosome, $position, $refBase, $altBase);
		my $sample = $tokenHash{'Sample'};
		$sample = $samplePatientHash{$sample} if(defined($samplePatientHash{$sample}));
		$geneGroupColumnValueHash{$tokenHash{'Gene name'}}->{$group}->{$tokenHash{'Mutation class'}}->{$sample}->{$variant} = 1 if($tokenHash{'Mutation class'} ne '');
	}
	close($reader);
}
my @groupList = sort keys %groupHash;
{
	my @tokenList = ('');
	foreach my $group (@groupList) {
		push(@tokenList, $group, ('') x (scalar(@columnList) - 1));
	}
	print join("\t", @tokenList), "\n";
}
{
	my @tokenList = ('gene');
	foreach my $group (@groupList) {
		push(@tokenList, @columnList);
	}
	print join("\t", @tokenList), "\n";
}
foreach my $gene (sort keys %geneGroupColumnValueHash) {
	my @tokenList = ($gene);
	foreach my $group (@groupList) {
		if(defined(my $columnValueHash = $geneGroupColumnValueHash{$gene}->{$group})) {
			foreach my $mutation (grep {defined($columnValueHash->{$_})} grep {$_ ne 'silent'} @mutationList) {
				foreach my $sample (keys %{$columnValueHash->{$mutation}}) {
					$columnValueHash->{'NS/SS/I'}->{$sample}->{$_} = 1 foreach(keys %{$columnValueHash->{$mutation}->{$sample}});
				}
			}
			foreach my $mutation (grep {defined($columnValueHash->{$_})} @mutationList, 'NS/SS/I') {
				foreach my $sample (keys %{$columnValueHash->{$mutation}}) {
					$columnValueHash->{"$mutation.variant"} += scalar(keys %{$columnValueHash->{$mutation}->{$sample}});
				}
			}
			push(@tokenList, map {defined($_) ? (ref($_) ? scalar(keys %$_) : $_) : ''} map {$columnValueHash->{$_}} @columnList);
		} else {
			push(@tokenList, ('') x scalar(@columnList));
		}
	}
	print join("\t", @tokenList), "\n";
}

sub getAltBaseAllele {
	my ($altBase) = @_;
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	return ($altBase, 1) if(scalar(@altBaseList) == 1);
	foreach my $altBaseIndex (0 .. $#altBaseList) {
		return ($1, $altBaseIndex + 1) if($altBaseList[$altBaseIndex] =~ /^\[(.*)\]$/);
	}
}
