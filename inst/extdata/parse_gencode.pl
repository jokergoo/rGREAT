use strict;

my $file = shift @ARGV;
open F, "gunzip -d $file -c |";

while(my $line = <F>) {
	my @line = split "\t", $line;

	if($line[2] ne "gene") {
		next;
	}

	if($line[8] =~/protein_coding/) {

		my @meta = split "; ", $line[8];
		my $gene_id = $meta[0];
		$gene_id =~s/gene_id//;
		$gene_id =~s/"//g;
		$gene_id =~s/ //g;

		print "$line[0]\t$line[3]\t$line[4]\t$line[6]\t$gene_id\n";
	}
}