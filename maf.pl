# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(all);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	't=s' => \(my $transcriptIdFile = ''),
	'd=s' => \(my $delimiter = ','),
	's=s' => \(my $sample = ''),
);
my %transcriptIdHash = ();
if($transcriptIdFile ne '') {
	open(my $reader, ($transcriptIdFile =~ /\.gz$/ ? "gzip -dc $transcriptIdFile |" : $transcriptIdFile));
	while(my $line = <$reader>) {
		chomp($line);
		$transcriptIdHash{$line} = 1;
	}
	close($reader);
}
my @mafColumnList = ('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'Annotation_Transcript', 'cDNA_Change', 'Protein_Change');
print join("\t", @mafColumnList), "\n";
my (@variantFileList) = @ARGV;
open(my $writer, "| sort -u");
foreach my $variantFile (@variantFileList) {
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
		next if($altBase eq '*');
		my %tokenListHash = ();
		@tokenListHash{@columnList} = map {[split(/$delimiter/, $_, -1)]} @tokenHash{@columnList};
		my $length = scalar(my @transcriptIdList = @{$tokenListHash{'Transcript ID'}});
		my @indexList = ();
		foreach my $index (0 .. $length - 1) {
			push(@indexList, $index) if($transcriptIdHash{$transcriptIdList[$index]});
		}
		unless(@indexList) {
			my @prefixNumberIndexList = ();
			foreach my $index (0 .. $length - 1) {
				(my $prefix = $transcriptIdList[$index]) =~ s/([0-9]+).*$//;
				push(@prefixNumberIndexList, [$prefix, $1, $index]);
			}
			@prefixNumberIndexList = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @prefixNumberIndexList;
			foreach(@prefixNumberIndexList) {
				if($_->[0] eq $prefixNumberIndexList[0]->[0] && $_->[1] == $prefixNumberIndexList[0]->[1]) {
					push(@indexList, $_->[2]);
				} else {
					last;
				}
			}
		}
		$tokenListHash{$_} = [($tokenHash{$_}) x $length] foreach(grep {scalar(@{$tokenListHash{$_}}) != $length} @columnList);
		foreach my $index (@indexList) {
			my %tokenHash = ();
			@tokenHash{@columnList} = map {$_->[$index]} @tokenListHash{@columnList};
			$tokenHash{'Hugo_Symbol'} = $tokenHash{'Gene name'};
			@tokenHash{'Start_Position', 'End_Position'} = ($position, $position + length($refBase) - 1);
			@tokenHash{'Reference_Allele', 'Tumor_Seq_Allele2'} = ($refBase, $altBase);
			if($tokenHash{'Mutation class'} eq 'silent') {
				$tokenHash{'Variant_Classification'} = 'Silent';
			} elsif($tokenHash{'Mutation class'} eq 'missense') {
				$tokenHash{'Variant_Classification'} = 'Missense_Mutation';
			} elsif($tokenHash{'Mutation class'} eq 'inframe' && $refBase =~ /$altBase/) {
				$tokenHash{'Variant_Classification'} = 'In_Frame_Del';
			} elsif($tokenHash{'Mutation class'} eq 'inframe' && $altBase =~ /$refBase/) {
				$tokenHash{'Variant_Classification'} = 'In_Frame_Ins';
			} elsif($tokenHash{'Mutation class'} eq 'frameshift' && $refBase =~ /$altBase/) {
				$tokenHash{'Variant_Classification'} = 'Frame_Shift_Del';
			} elsif($tokenHash{'Mutation class'} eq 'frameshift' && $altBase =~ /$refBase/) {
				$tokenHash{'Variant_Classification'} = 'Frame_Shift_Ins';
			} elsif($tokenHash{'Mutation class'} eq 'nonsense') {
				$tokenHash{'Variant_Classification'} = 'Nonsense_Mutation';
			} elsif($tokenHash{'Mutation class'} eq 'readthrough') {
				$tokenHash{'Variant_Classification'} = 'Nonstop_Mutation';
			} elsif($tokenHash{'Mutation class'} eq 'startcodon') {
				$tokenHash{'Variant_Classification'} = 'Translation_Start_Site';
			} elsif($tokenHash{'Mutation class'} eq 'splicing') {
				$tokenHash{'Variant_Classification'} = 'Splice_Site';
			} elsif($tokenHash{'Region type'} eq 'utr5') {
				$tokenHash{'Variant_Classification'} = "5'UTR";
			} elsif($tokenHash{'Region type'} eq 'utr3') {
				$tokenHash{'Variant_Classification'} = "3'UTR";
			} elsif($tokenHash{'Region type'} eq 'intron') {
				$tokenHash{'Variant_Classification'} = 'Intron';
			} elsif($tokenHash{'Region type'} eq 'non_coding_exon') {
				$tokenHash{'Variant_Classification'} = 'RNA';
			} elsif($tokenHash{'Region type'} eq 'non_coding_intron') {
				$tokenHash{'Variant_Classification'} = 'Intron';
			}
			$tokenHash{'Variant_Type'} = getVariantType($refBase, $altBase);
			$tokenHash{'Tumor_Sample_Barcode'} = $sample ne '' ? $sample : $tokenHash{'Sample'};
			$tokenHash{'Annotation_Transcript'} = $tokenHash{'Transcript ID'};
			$tokenHash{'cDNA_Change'} = $tokenHash{'Transcript variation nomenclature'};
			$tokenHash{'Protein_Change'} = $tokenHash{'Protein variation nomenclature'};
			if(all {defined($_)} @tokenHash{@mafColumnList}) {
				print $writer join("\t", @tokenHash{@mafColumnList}), "\n";
			} else {
#				print STDERR join("\t", map {defined($_) ? $_ : ''} @tokenHash{@mafColumnList}), "\n";
			}
		}
		unless(@indexList) {
			$tokenHash{'Hugo_Symbol'} = 'Unknown';
			@tokenHash{'Start_Position', 'End_Position'} = ($position, $position + length($refBase) - 1);
			@tokenHash{'Reference_Allele', 'Tumor_Seq_Allele2'} = ($refBase, $altBase);
			$tokenHash{'Variant_Classification'} = 'IGR';
			$tokenHash{'Variant_Type'} = getVariantType($refBase, $altBase);
			$tokenHash{'Tumor_Sample_Barcode'} = $sample ne '' ? $sample : $tokenHash{'Sample'};
			$tokenHash{'Annotation_Transcript'} = '';
			$tokenHash{'cDNA_Change'} = '';
			$tokenHash{'Protein_Change'} = '';
			if(all {defined($_)} @tokenHash{@mafColumnList}) {
				print $writer join("\t", @tokenHash{@mafColumnList}), "\n";
			} else {
#				print STDERR join("\t", map {defined($_) ? $_ : ''} @tokenHash{@mafColumnList}), "\n";
			}
		}
	}
	close($reader);
}
close($writer);

sub getAltBaseAllele {
	my ($altBase) = @_;
	my @altBaseList = $altBase eq '' ? ('') : split(/,/, $altBase, -1);
	return ($altBase, 1) if(scalar(@altBaseList) == 1);
	foreach my $altBaseIndex (0 .. $#altBaseList) {
		return ($1, $altBaseIndex + 1) if($altBaseList[$altBaseIndex] =~ /^\[(.*)\]$/);
	}
}

sub getVariantType {
	my ($refBase, $altBase) = @_;
	if($refBase =~ /^[ACGT]*$/ && $altBase =~ /^[ACGT]*$/ && $refBase ne $altBase) {
		return 'SNP' if(length($refBase) == 1 && length($altBase) == 1);
		return 'DNP' if(length($refBase) == 2 && length($altBase) == 2);
		return 'TNP' if(length($refBase) == 3 && length($altBase) == 3);
		return 'ONP' if(length($refBase) == length($altBase));
		return 'INS' if($altBase =~ /$refBase/);
		return 'DEL' if($refBase =~ /$altBase/);
	}
}
