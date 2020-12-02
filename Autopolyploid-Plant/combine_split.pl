open OUT,">R2-$ARGV[1].txt" or die "can't open the output file";
foreach $i(1..$ARGV[2]){
	$num=0;
	open IN,"R2-split$i-$ARGV[1].txt" or die "can't open the input file";
	while (<IN>){
		chomp();
		print OUT "$_\n";
	}
}