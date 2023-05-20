# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use List::Util qw(all sum min max reduce);
use Getopt::Long qw(:config no_ignore_case);

chomp(my $hostname = `hostname`);
my $temporaryDirectory = $ENV{'TMPDIR'};
$temporaryDirectory = '/tmp' unless($temporaryDirectory);
GetOptions(
	'h' => \(my $help = ''),
	't=s' => \$temporaryDirectory,
	'p=i' => \(my $threads = 1),
	'numberPerThread=i' => \(my $numberPerThread = 10000),
	'c=s' => \(my $columnNameFile = ''),
	'N=s' => \(my $normalSample = ''),
	'T=s' => \(my $tumorSample = ''),
	'P' => \(my $passFilter = ''),
	'I' => \(my $doNotParseInfoField = ''),
	'G' => \(my $doNotParseGenotypeField = ''),
	'S' => \(my $doNotPrintSingleGenotypeVariant = ''),
	'Z' => \(my $doNotPrintHomozygousReferenceGenotype = ''),
	'V' => \(my $doNotPrintUnmatchedVariant = ''),
	'C' => \(my $doNotPrintCommonVariant = ''),
	'U' => \(my $doNotUnphaseGenotype = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vcf.table.pl [options] input.vcf column=name [...] > variant.txt

Options: -h       display this help message
         -t DIR   temporary directory [$temporaryDirectory]
         -p INT   threads [$threads]
         -c FILE  column name file
         -N STR   normal sample
         -T STR   tumor sample
         -P       pass filter
         -I       do not parse INFO field
         -G       do not parse Genotype field
         -S       do not print single genotype variant
         -Z       do not print homozygous reference genotype
         -V       do not print unmatched variant
         -C       do not print common variant
         -U       do not unphase genotype

EOF
}
{
	my $parentPid = $$;
	my %pidHash = ();
	my $writer;
	my $parentWriter;
	sub forkPrintParentWriter {
		($parentWriter) = @_;
	}
	sub forkPrintSubroutine {
		my ($subroutine, @arguments) = @_;
		if(my $pid = fork()) {
			$pidHash{$pid} = 1;
		} else {
			open($writer, "> $temporaryDirectory/fork.$hostname.$parentPid.$$");
			$subroutine->(@arguments);
			close($writer);
			exit(0);
		}
		forkPrintWait($threads);
	}
	sub forkPrintWait {
		my ($number) = (@_, 1);
		while(scalar(keys %pidHash) >= $number) {
			my $pid = wait();
			if($pidHash{$pid}) {
				open(my $reader, "$temporaryDirectory/fork.$hostname.$parentPid.$pid");
				if(defined($parentWriter)) {
					print $parentWriter $_ while(<$reader>);
				} else {
					print $_ while(<$reader>);
				}
				close($reader);
				system("rm $temporaryDirectory/fork.$hostname.$parentPid.$pid");
				delete $pidHash{$pid};
			}
		}
	}
	sub forkPrint {
		if(defined($writer)) {
			print $writer @_;
		} elsif(defined($parentWriter)) {
			print $parentWriter @_;
		} else {
			print @_;
		}
	}
}
my ($vcfFile, @columnNameList) = @ARGV;
if(@columnNameList) {
	@columnNameList = map {[split(/=/, $_, 2)]} @columnNameList;
} elsif($columnNameFile ne '') {
	open(my $reader, $columnNameFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my @tokenList = split(/\t/, $line, -1);
		push(@columnNameList, [$tokenList[0], @tokenList[1 .. $#tokenList]]);
	}
	close($reader);
}
my @columnList = map {$_->[0]} @columnNameList;
my @sampleColumnList = ();
my @tumorSampleColumnList = ();
my @normalSampleColumnList = ();
my @functionColumnList = ();
my @titleList = ();
foreach my $column (@columnList) {
	if($column =~ /^SAMPLE\.(.*)$/) {
		push(@sampleColumnList, $1);
	} elsif($column =~ /^tumor\.(.*)$/) {
		push(@tumorSampleColumnList, $1);
	} elsif($column =~ /^normal\.(.*)$/) {
		push(@normalSampleColumnList, $1);
	} elsif($column =~ /^(.*)\((.*)\)$/) {
		push(@functionColumnList, $column);
	} elsif($column =~ /^([^.]*)\.(.*)$/) {
		push(@titleList, $1);
	}
}
@titleList = unique(@titleList) if(@titleList);
my @sampleList = ();
my %sampleIndexHash = ();
{
	open(my $reader, ($vcfFile =~ /\.gz$/ ? "gzip -dc $vcfFile |" : $vcfFile));
	my @lineList = ();
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO//) {
			if($line =~ s/^\tFORMAT\t//) {
				@sampleList = split(/\t/, $line, -1);
				@sampleIndexHash{@sampleList} = 0 .. $#sampleList;
				if($doNotParseGenotypeField eq '') {
					foreach my $index (0 .. $#columnNameList) {
						if($columnNameList[$index]->[0] =~ /^SAMPLE\.(.*)$/) {
							if(defined($columnNameList[$index]->[1])) {
								$columnNameList[$index]->[1] = join("\t", map {"$_ $columnNameList[$index]->[1]"} @sampleList);
							} else {
								$columnNameList[$index]->[1] = join("\t", map {"$_.$1"} @sampleList);
							}
						}
					}
				}
			}
			print '#', join("\t", map {$_->[-1]} @columnNameList), "\n";
			next;
		}
		if($line =~ /^#/) {
			$normalSample = $1 if($normalSample eq '' && $line =~ / --normal-sample (\S+) /);
			$tumorSample = $1 if($tumorSample eq '' && $line =~ / --tumor-sample (\S+) /);
			next;
		}
		push(@lineList, $line);
		if(scalar(@lineList) >= $numberPerThread) {
			if($threads == 1) {
				printTable(@lineList);
			} else {
				forkPrintSubroutine(\&printTable, @lineList);
			}
			@lineList = ();
		}
	}
	if(@lineList) {
		if($threads == 1) {
			printTable(@lineList);
		} else {
			forkPrintSubroutine(\&printTable, @lineList);
		}
		@lineList = ();
	}
	forkPrintWait();
	close($reader);
}

sub printTable {
	foreach my $line (@_) {
		my %tokenHash = ();
		(@tokenHash{'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}, my @genotypeList) = map {$_ eq '.' ? '' : $_} split(/\t/, $line, -1);
		next if($passFilter eq '' && $tokenHash{'FILTER'} ne 'PASS');
		$tokenHash{'END'} = $tokenHash{'POS'} + length($tokenHash{'REF'}) - 1;
		my $altBase = $tokenHash{'ALT'};
		my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
		my @titleTokenHashListHashList = ();
		if($doNotParseInfoField eq '') {
			my %titleTokenHashListHash = ();
			foreach my $columnValue (split(';', $tokenHash{'INFO'})) {
				if($columnValue =~ /^([^=]*)=(.*)$/) {
					my ($column, $value) = ($1, $2);
					if($column =~ /^([^.]*)\.(.*)$/) {
						my ($title, $index) = ($1, 0);
						($column, $index) = ($1, $2 - 1) if($column =~ /^(.*)_([0-9]+)$/);
						$titleTokenHashListHash{$title}->[$index]->{$column} = $value;
					} else {
						$tokenHash{$column} = $value;
					}
				} else {
					$tokenHash{$columnValue} = $columnValue unless(defined($tokenHash{$columnValue}));
				}
			}
			if(scalar(@altBaseList) > 1) {
				foreach my $altBaseIndex (0 .. $#altBaseList) {
					foreach my $title (keys %titleTokenHashListHash) {
						$titleTokenHashListHashList[$altBaseIndex]->{$title} = [grep {!defined($_->{"$title.allele"}) || $_->{"$title.allele"} == $altBaseIndex + 1} @{$titleTokenHashListHash{$title}}];
					}
				}
			} else {
				push(@titleTokenHashListHashList, \%titleTokenHashListHash);
			}
		}
		my @tokenHashList = ();
		if(defined($tokenHash{'FORMAT'}) && $doNotParseGenotypeField eq '') {
			my @formatList = split(/:/, $tokenHash{'FORMAT'}, -1);
			foreach my $index (0 .. $#sampleList) {
				$tokenHashList[$index]->{'SAMPLE'} = $sampleList[$index];
				@{$tokenHashList[$index]}{@formatList} = split(/:/, $genotypeList[$index], -1);
			}
		}
		if($doNotPrintSingleGenotypeVariant) {
			my %genotypeHash = ();
			$genotypeHash{$_} = 1 foreach(map {unphaseGenotype($_->{'GT'})} @tokenHashList);
			next if(scalar(keys %genotypeHash) == 1);
		}
		if($doNotPrintHomozygousReferenceGenotype) {
			@tokenHashList = grep {!isHomozygousReferenceGenotype($_->{'GT'})} @tokenHashList;
			next if(scalar(@tokenHashList) == 0);
		}
		if($doNotUnphaseGenotype eq '') {
			$_->{'GT'} = unphaseGenotype($_->{'GT'}) foreach(@tokenHashList);
		}
		foreach my $column (@sampleColumnList) {
			my %sampleValueHash = map {$_->{'SAMPLE'} => $_->{$column}} @tokenHashList;
			$tokenHash{"SAMPLE.$column"} = join("\t", map {defined($_) ? $_ : ''} @sampleValueHash{@sampleList});
		}
		foreach my $altBaseIndex (0 .. $#altBaseList) {
			my %titleTokenHashListHash = ();
			if(@titleTokenHashListHashList) {
				%titleTokenHashListHash = %{$titleTokenHashListHashList[$altBaseIndex]};
				foreach my $column (@functionColumnList) {
					if($column =~ /^(.*)\((.*)\)$/) {
						my ($function, $column1) = ($1, $2);
						if($column1 =~ /^(.*)\.(.*)=(.*)\.(.*);(.*)$/) {
							my ($title1, $column1, $title2, $column2, $column3) = ($1, "$1.$2", $3, "$3.$4", "$3.$5");
							if(defined(my $tokenHashList1 = $titleTokenHashListHash{$title1})) {
								if(defined(my $tokenHashList2 = $titleTokenHashListHash{$title2})) {
									my %tokenHash = ();
									foreach(grep {defined($_->{$column2}) && defined($_->{$column3})} @$tokenHashList2) {
										push(@{$tokenHash{$_->{$column2}}}, $_->{$column3});
									}
									foreach(keys %tokenHash) {
										$tokenHash{$_} = executeFunction($function, @{$tokenHash{$_}});
									}
									foreach(grep {defined($_->{$column1})} @$tokenHashList1) {
										$_->{$column} = $tokenHash{$_->{$column1}};
									}
								}
							}
						} elsif($column1 =~ /^(.*)\.(.*)$/) {
							my ($title1, $column1) = ($1, "$1.$2");
							if(defined(my $tokenHashList1 = $titleTokenHashListHash{$title1})) {
								$tokenHash{$column} = executeFunction($function, grep {defined} map {$_->{$column1}} @$tokenHashList1);
							} else {
								$tokenHash{$column} = '';
							}
						}
					}
				}
			}
			if(scalar(@altBaseList) > 1) {
				$tokenHash{'ALT'} = join(',', map {$_ == $altBaseIndex ? "[$altBaseList[$_]]" : $altBaseList[$_]} 0 .. $#altBaseList);
			}
			foreach my $tokenHash (@tokenHashList) {
				if($tokenHash->{'AD'}) {
					my @alleleDepthList = split(/,/, $tokenHash->{'AD'});
					$tokenHash->{'refAD'} = $alleleDepthList[0];
					$tokenHash->{'altAD'} = $alleleDepthList[$altBaseIndex + 1];
					if(defined(my $depth = $tokenHash->{'DP'})) {
						$tokenHash->{'AF'} = $depth > 0 ? $tokenHash->{'altAD'} / $depth : "$tokenHash->{'altAD'}/$depth";
					}
				}
			}
			if($normalSample ne '' && defined($sampleIndexHash{$normalSample})) {
				foreach my $column (@normalSampleColumnList) {
					$tokenHash{"normal.$column"} = $tokenHashList[$sampleIndexHash{$normalSample}]->{$column};
				}
			}
			if($tumorSample ne '' && defined($sampleIndexHash{$tumorSample})) {
				foreach my $column (@tumorSampleColumnList) {
					$tokenHash{"tumor.$column"} = $tokenHashList[$sampleIndexHash{$tumorSample}]->{$column};
				}
			}
			my $pid = open2(my $reader, my $writer, 'sort -u');
			if($doNotPrintUnmatchedVariant || $doNotPrintCommonVariant) {
				if(my @matchedTokenHashList = grep {grep {$_ == $altBaseIndex + 1} split(/[\/|]/, $_->{'GT'})} @tokenHashList) {
					unless($doNotPrintCommonVariant && scalar(@matchedTokenHashList) == scalar(@tokenHashList)) {
						printLines($writer, \%titleTokenHashListHash, \@matchedTokenHashList, \%tokenHash, 0);
					}
				}
			} else {
				printLines($writer, \%titleTokenHashListHash, \@tokenHashList, \%tokenHash, 0);
			}
			close($writer);
			forkPrint($_) while(<$reader>);
			close($reader);
			waitpid($pid, 0);
		}
	}
}

sub printLines {
	my ($writer, $titleTokenHashListHash, $tokenHashList, $tokenHash, $titleIndex) = @_;
	if($titleIndex < scalar(@titleList)) {
		if(defined($_ = $titleTokenHashListHash->{$titleList[$titleIndex]}) && (my @tokenHashList = @$_)) {
			printLines($writer, $titleTokenHashListHash, $tokenHashList, {%$tokenHash, %$_}, $titleIndex + 1) foreach(@tokenHashList);
		} else {
			printLines($writer, $titleTokenHashListHash, $tokenHashList, $tokenHash, $titleIndex + 1);
		}
	} else {
		if(@$tokenHashList) {
			print $writer join("\t", map {defined($_) ? $_ : ''} @$_{@columnList}), "\n" foreach(map {{%$tokenHash, %$_}} @$tokenHashList);
		} else {
			print $writer join("\t", map {defined($_) ? $_ : ''} @$tokenHash{@columnList}), "\n";
		}
	}
}

sub executeFunction {
	my ($function, @tokenList) = @_;
	if(@tokenList) {
		if($function eq '') {
			if(all {$_ eq $tokenList[0]} @tokenList[1 .. $#tokenList]) {
				return $tokenList[0];
			} else {
				return join(',', @tokenList);
			}
		}
		return join(',', unique(@tokenList)) if($function eq 'unique');
		return sum(@tokenList) if($function eq 'sum');
		return min(@tokenList) if($function eq 'min');
		return max(@tokenList) if($function eq 'max');
		return reduce {abs($a) < abs($b) ? $a : $b} @tokenList if($function eq 'minabs');
		return reduce {abs($a) > abs($b) ? $a : $b} @tokenList if($function eq 'maxabs');
	} else {
		return '';
	}
}

sub unphaseGenotype {
	return join('/', sort {$a <=> $b} grep {$_ ne '.'} split(/[\/|]/, $_[0])) if(scalar(@_) == 1);
	return map {unphaseGenotype($_)} @_;
}

sub isHomozygousReferenceGenotype {
	my ($genotype) = @_;
	my @alleleList = unique(grep {$_ ne '.'} split(/[\/|]/, $genotype));
	return scalar(@alleleList) == 1 && $alleleList[0] eq '0';
}

sub unique {
	my @sortedTokenList = sort @_;
	@sortedTokenList = @sortedTokenList[0, grep {$sortedTokenList[$_ - 1] ne $sortedTokenList[$_]} 1 .. $#sortedTokenList] if(@sortedTokenList);
	return @sortedTokenList;
}
