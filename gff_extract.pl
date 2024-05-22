#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(all);
use Getopt::Long qw(:config no_ignore_case);

my @filterList = ();
GetOptions(
	'h' => \(my $help = ''),
	'f=s' => \@filterList,
	'b=s' => \(my $booleanOperator = 'and'),
	'E' => \(my $noEmpty = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl gff_extract.pl [options] gene.gff [column ...] > extract.txt

Options: -h       display this help message
         -f STR   filter e.g. feature=exon
         -b STR   boolean operator, "and" or "or" [$booleanOperator]
         -E       no empty

EOF
}
@filterList = map {[split(/=/, $_, 2)]} @filterList;
my ($gffFile, @columnList) = @ARGV;
open(my $reader, ($gffFile =~ /\.gz$/) ? "gzip -dc $gffFile | grep -v '^#' | LC_ALL=C sort -t '\t' -k1,1 -k4,4n -k5,5n |" : "grep -v '^#' $gffFile | LC_ALL=C sort -t '\t' -k1,1 -k4,4n -k5,5n |");
my %tokenHashHash = ();
while(my $line = <$reader>) {
	chomp($line);
	my %tokenHash = ();
	@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
	my %attributeHash = ();
	$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^;= ]+)=([^;]+)(;|$)/g);
	$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^;" ]+) +"([^;"]+)"(;|$)/g);
	$tokenHash{'attribute'} = \%attributeHash;
	my @deleteIdList = ();
	foreach my $id (keys %tokenHashHash) {
		next if($tokenHashHash{$id}->{'chromosome'} eq $tokenHash{'chromosome'} && $tokenHash{'start'} <= $tokenHashHash{$id}->{'end'});
		next if(getAncestor($id) eq '');
		extract($id);
		push(@deleteIdList, $id);
	}
	delete $tokenHashHash{$_} foreach(@deleteIdList);
	if(defined(my $id = $tokenHash{'attribute'}->{'ID'})) {
		extract($id) if(defined($tokenHashHash{$id}));
		$tokenHashHash{$id} = \%tokenHash;
	} else {
		$tokenHashHash{$tokenHash{'attribute'}} = \%tokenHash;
	}
}
{
	my @deleteIdList = ();
	foreach my $id (keys %tokenHashHash) {
		extract($id);
		push(@deleteIdList, $id);
	}
	delete $tokenHashHash{$_} foreach(@deleteIdList);
}
close($reader);

sub extract {
	my ($id) = @_;
	if($booleanOperator eq 'and') {
		if(all {getColumnValue($id, $_->[0]) eq $_->[1]} @filterList) {
			my @valueList = map {getColumnValue($id, $_)} @columnList;
			return if($noEmpty && grep {$_ eq ''} @valueList);
			print join("\t", @valueList), "\n";
		}
	}
	if($booleanOperator eq 'or') {
		if(grep {getColumnValue($id, $_->[0]) eq $_->[1]} @filterList) {
			my @valueList = map {getColumnValue($id, $_)} @columnList;
			return if($noEmpty && grep {$_ eq ''} @valueList);
			print join("\t", @valueList), "\n";
		}
	}
}

sub getAncestor {
	my ($id) = @_;
	if(defined($tokenHashHash{$id})) {
		if(defined(my $parentId = $tokenHashHash{$id}->{'attribute'}->{'Parent'})) {
			return getAncestor($parentId);
		} else {
			return $id;
		}
	} else {
		return '';
	}
}

sub getColumnValue {
	my ($id, $column) = @_;
	foreach my $column (split(/,/, $column)) {
		if(defined(my $value = $tokenHashHash{$id}->{$column})) {
			return $value;
		}
		if(defined(my $value = getAttributeValue($id, $column))) {
			return $value;
		}
	}
	return '';
}

sub getAttributeValue {
	my ($id, $attribute) = @_;
	if(defined(my $value = $tokenHashHash{$id}->{'attribute'}->{$attribute})) {
		return $value;
	}
	if(defined(my $parentId = $tokenHashHash{$id}->{'attribute'}->{'Parent'})) {
		return getAttributeValue($parentId, $attribute);
	}
	return undef;
}
