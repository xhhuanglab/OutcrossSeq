$chr =$ARGV[0];
if($ARGV[1]%$ARGV[2]>=0.6*$ARGV[2]){
	$file_num=int($ARGV[1]/$ARGV[2])+1;
}else{
	$file_num=int($ARGV[1]/$ARGV[2]);
}
 $num=1;
open IN,"$ARGV[3].geno" or die "can't open $ARGV[3]";
while (<IN>){
	chomp();
	if ($.==1){
		$head=$_;
	}else{
		$hang++;
		$line=$_;
		if($hang==1 and $num<=$file_num){
			open OUT,">split$num-$ARGV[3].geno";
			print OUT "$head\n$line\n";
		}elsif($hang==$ARGV[2] and $num<$file_num){
			$hang=0;$num++;
			print OUT "$line\n";
		}else{
			print OUT "$line\n";
		}
	}
}
close IN;
close OUT;
