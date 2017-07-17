# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

my ($gbFile, $directory) = @ARGV;
system("mkdir -p $directory");
open(my $reader, $gbFile);
while(my $line = <$reader>) {
	if($line =~ /^LOCUS +([^ ]+) +/) {
		my $locusName = $1;
		open(my $writer, "> $directory/$locusName.gb");
		print $writer $line;
		while($line && $line ne "//\n") {
			print $writer ($line = <$reader>);
		}
		close($writer);
	}
}
close($reader);
