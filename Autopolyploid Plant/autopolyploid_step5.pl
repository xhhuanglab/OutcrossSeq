open OUT,">genotype.txt" ;
foreach $chr($ARGV[0]..$ARGV[1]){
	$in="Chr$chr-Imputation_marker_sort.txt";
	open IN,"$in" or die "can't open $in";
	while (<IN>){
		chomp();
		$line =$_;
		if ($.==1){
			if ($chr == $ARGV[0]){
				print OUT "$line\n";
			}
		}else{
			if($line=~/Chr(\d+)/){
				if($chr==$1+1-1){
					print OUT "$line\n";
				}
			}
		}
	}
	close IN;
}
close OUT;