open IN1,"sample.list";
while (<IN1>){
	chomp();
	push @sample,$_;
}
close IN1;
open IN,"$ARGV[0]" or die "can't open the input file";
while (<IN>){ 
	chomp();
	$line =$_;
	if ($line=~/^#CHROM/){
		@line = split/\s+/,$line;
		foreach $j(@sample){
			foreach $i(8..$#line){
				if ($j eq $line[$i]){
					push @i,$i;
				}
			}
		}
	}elsif ($line=~/Chr(\d+)\s+(\d+)\s+([A T G C])\s([A T G C])/){
		$key="$1	$2	$3	$4";
		$reads_num{$key}=$line;
	}
} 
open IN2,"$ARGV[1].txt" or die "can't open the input file";
open OUT,">$ARGV[1]-readsNum.txt";
while (<IN2>){
	$array=$_;
	if ($array=~/Chr(\d+)\s+(\d+)\s+([A T G C])\s([A T G C])/){
		chomp();
		$key="$1	$2	$3	$4";
		@marker=split/\s+/,$reads_num{$key};
		print OUT "$marker[0]	$marker[1]	$marker[2]	$marker[3]	$marker[4]	$marker[5]	$marker[6]	$marker[7]";
		foreach (@i){
			print OUT "\t$marker[$_]";
		}
		print OUT "\n";
	}elsif($array=~/^\n/){
		print OUT "\n";
	}else{
		print OUT "none\n";
	}
}
close IN;
close IN1;
close IN2;
close OUT;
