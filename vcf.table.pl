# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum min max reduce);
use Getopt::Long;

GetOptions('p' => \(my $pass = 0));
my ($vcfFile, @columnList) = @ARGV;
@columnList = map {split(/,/, $_)} @columnList;
my @titleList = sort map {$_->[0]} grep {scalar(@$_) > 1} map {[split(/\./, $_)]} grep {$_ !~ /^(.*)\((.*)\)$/} @columnList;
@titleList = @titleList[0, grep {$titleList[$_ - 1] ne $titleList[$_]} 1 .. $#titleList] if(@titleList);
my $headerLine = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
my @sampleList = ();
open(my $reader, $vcfFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ s/^#$headerLine\t// && (@sampleList = split(/\t/, $line)));
	next if($line =~ /^#/);
	my %tokenHash = ();
	(@tokenHash{split(/\t/, $headerLine)}, my @genotypeList) = map {$_ eq '.' ? '' : $_} split(/\t/, $line);
	next if($pass == 0 && $tokenHash{'FILTER'} ne 'PASS');
	my @tokenHashList = ();
	for(my $index = 0; $index < scalar(@sampleList); $index++) {
		$tokenHashList[$index] = {'SAMPLE' => $sampleList[$index]};
		@{$tokenHashList[$index]}{split(':', $tokenHash{'FORMAT'})} = split(':', $genotypeList[$index]);
	}
	next if($pass == 0 && grep {$_->{'GT'} eq './.'} @tokenHashList);
	my %tokenHashListHash = ();
	foreach my $keyValue (split(';', $tokenHash{'INFO'})) {
		if($keyValue =~ /^([^=]*)=(.*)$/) {
			my ($key, $value) = ($1, $2);
			if($key =~ /^([^.]*)\.(.*)$/) {
				my ($title, $index) = ($1, 0);
				($key, $index) = ($1, $2 - 1) if($key =~ /^(.*)_([0-9]+)$/);
				$tokenHashListHash{$title}->[$index]->{$key} = $value;
			} else {
				$tokenHash{$key} = $value;
			}
		} else {
			$tokenHash{$keyValue} = $keyValue;
		}
	}
	foreach my $column (@columnList) {
		if($column =~ /^(.*)\((.*)\)$/) {
			my ($function, $column1) = ($1, $2);
			if($column1 =~ /^(.*)\.(.*)=(.*)\.(.*);(.*)$/) {
				my ($title1, $column1, $title2, $column2, $column3) = ($1, "$1.$2", $3, "$3.$4", "$3.$5");
				next unless(defined($tokenHashListHash{$title1}));
				next unless(defined($tokenHashListHash{$title2}));
				my %valueHash = ();
				foreach(grep {defined($_->{$column2})} @{$tokenHashListHash{$title2}}) {
					push(@{$valueHash{$_->{$column2}}}, $_->{$column3});
				}
				foreach(keys %valueHash) {
					$valueHash{$_} = executeFunction($function, grep {defined} @{$valueHash{$_}});
				}
				foreach(grep {defined($_->{$column1})} @{$tokenHashListHash{$title1}}) {
					$_->{$column} = $valueHash{$_->{$column1}};
				}
			} elsif($column1 =~ /^(.*)\.(.*)$/) {
				my ($title1, $column1) = ($1, "$1.$2");
				next unless(defined($tokenHashListHash{$title1}));
				$tokenHash{$column} = executeFunction($function, grep {defined} map {$_->{$column1}} @{$tokenHashListHash{$title1}});
			}
		}
	}
	open(my $writer, "| sort -u");
	printLines($writer, \%tokenHashListHash, \@tokenHashList, \%tokenHash, 0);
	close($writer);
}
close($reader);

sub printLines {
	my ($writer, $tokenHashListHash, $tokenHashList, $tokenHash, $index) = @_;
	if($index < scalar(@titleList)) {
		if(defined($tokenHashListHash->{$titleList[$index]})) {
			printLines($writer, $tokenHashListHash, $tokenHashList, {%$tokenHash, %$_}, $index + 1) foreach(@{$tokenHashListHash->{$titleList[$index]}});
		} else {
			printLines($writer, $tokenHashListHash, $tokenHashList, $tokenHash, $index + 1);
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
	my ($function, @valueList) = @_;
	return undef if(scalar(@valueList) == 0);
	return join(',', @valueList) if($function eq '');
	return join(',', unique(@valueList)) if($function eq 'unique');
	return sum(@valueList) if($function eq 'sum');
	return min(@valueList) if($function eq 'min');
	return max(@valueList) if($function eq 'max');
	return reduce {abs($a) < abs($b) ? $a : $b} @valueList if($function eq 'minabs');
	return reduce {abs($a) > abs($b) ? $a : $b} @valueList if($function eq 'maxabs');
}

sub unique {
	my @sorted = sort @_;
	return @sorted[0, grep {$sorted[$_ - 1] ne $sorted[$_]} 1 .. $#sorted];
}
