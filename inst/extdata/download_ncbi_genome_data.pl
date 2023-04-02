use strict;

while(my $line = <>) {
	chomp $line;
	my @lines = split "\t", $line;
	next if($lines[2] ne "gene");
	next if($lines[8] !~/protein_coding/);
	
	$lines[8] =~s/^.*Dbxref=GeneID:(\d+).*$/$1/;

	print "$lines[0]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[8]\n";
}
