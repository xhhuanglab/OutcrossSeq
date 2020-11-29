my %windows_num;my ($chr,$filename,@file,$marker);
$windows_num{0}=0;
foreach $chr($ARGV[0]..$ARGV[1]){
	$filename = "$chr".'_1modify.parent1';
    open IN1,"$filename" or die "can't open the input file 1";
    @file = <IN1>;
    $windows_num{$chr} = $#file+1+$windows_num{$chr-1};
    $marker = $#file+1;
}
my ($pos_star,$mark_num,$pos);
foreach $chr($ARGV[0]..$ARGV[1]){
	$pos_star = $ARGV[3]/2000000;
	foreach $mark_num($windows_num{$chr-1}+1..$windows_num{$chr}){
		$pos = $pos_star+($ARGV[3]/1000000)*($mark_num-$windows_num{$chr-1}-1);
		$chr[$mark_num]=$chr;
		$pos[$mark_num]=$pos;
	}
}

my ($chromsome,$sample,$filename1,$n,@line,@geno,$filename2,@array,$abscissa);
foreach $chromsome($ARGV[0]..$ARGV[1]){
	foreach $sample(1..$ARGV[2]){
		$filename1 = "$chromsome".'_'."$sample".'modify.parent1';
		open IN2,"$filename1" or die "can't open the input file 1";
		while (<IN2>){
			chomp();
			@line = split/\s+/,$_;
			$abscissa = $line[0]+$windows_num{$chromsome-1};
			$geno[$abscissa][$sample*2-1] = $line[1];
		}
	}
}
$n = 0;
foreach $chromsome($ARGV[0]..$ARGV[1]){
	foreach $sample(1..$ARGV[2]){
		$filename2 = "$chromsome".'_'."$sample".'modify.parent2';
		open IN3,"$filename2" or die "can't open the input file 2";
		$n=0;
		while (<IN3>){
			chomp(); 
			$n++;
			@array = split/\s+/,$_;
			$abscissa = $array[0]+$windows_num{$chromsome-1};
			$geno[$abscissa][$sample*2] = $array[1];
		}


	}
}
open OUT,">genotype.txt" or die "can't open the output file";
print OUT "Chr	Pos";
open IN4,"sample.list" or die "can't open the input file 4";
while (<IN4>){
	chomp();
	print OUT "\t$_";
}
print OUT "\n";	
my ($i,$j);my $sample_num = $ARGV[2]*2;
foreach $i(1..$windows_num{$ARGV[1]}){
	print OUT "$chr[$i]\t$pos[$i]";
    foreach  $j(1..$sample_num){
    	if ($geno[$i][$j] eq "A"){
            print OUT "\tA";
        }elsif($geno[$i][$j] eq "B"){
        	print OUT "\tB";
        }elsif($geno[$i][$j] eq "C"){
        	print OUT "C";
	}elsif($geno[$i][$j] eq "D"){
        	print OUT "D";
        }
     }
     print OUT "\n";
}
close IN1;
close IN2;
close IN3;
close IN4;
close OUT;