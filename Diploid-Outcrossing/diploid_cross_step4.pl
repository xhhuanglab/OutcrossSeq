$chr = $ARGV[1];
$hang=0;$chr_hang[0]=0;
foreach $chromosome($ARGV[0]..$chr){
	$chr_hang=0;
	$file="Imputation-Chr$chromosome.genotype";
	open IN,"$file" or die "can't open the in file";
	while (<IN>){
		chomp();
		@line = split/\s+/,$_;
		$sample_num=$#line-1;
		$hang++;$chr_hang++;
		$marker[$hang][0]=$chromosome;
		foreach (1..$#line){
			$marker[$hang][$_]=$line[$_];
		}
	}
	$chr_hang[$chromosome]=$chr_hang+$chr_hang[$chromosome-1];
	close IN;
}

open IN1,"sample.list" or die "can't open sample.list";
open OUT,">$ARGV[2].geno";
print OUT "Chr	Pos";
while (<IN1>){
	chomp();
	print OUT "\t$_";
}
print OUT "\n";
$chromosome=$ARGV[0];
$marker[$0][1]=0;
foreach $i(1..$hang){
	if($i<=$chr_hang[$chromosome]){
		$pos=(($marker[$i][1]-$marker[$i-1][1])/2+$marker[$i-1][1]+1)/1000000;$pos=sprintf "%.3f",$pos;
		print OUT "$chromosome	$pos";
		foreach (2..$#line){
			print OUT "\t$marker[$i][$_]";
		}
		print OUT "\n";
	}else{
		$chromosome++;$pos=($marker[$i][1]+1)/2000000;$pos=sprintf "%.3f",$pos;
		print OUT "$chromosome	$pos";
		foreach (2..$#line){
			print OUT "\t$marker[$i][$_]";
		}
		print OUT "\n";
	}
}
close IN1;
close OUT;
