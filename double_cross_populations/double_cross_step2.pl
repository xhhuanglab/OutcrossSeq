use List::Util qw/max min/;
my ($line,@line,$hang_num,$chr,@chr,$key,$pos);our (%s4,%s5,%s6,%s7,%chr_star,%chr_end);
open IN1,"$ARGV[0]" or die "can't open the input file 2";
while (<IN1>){
	chomp();
	$line = $_;
	if ($.>1){
		$hang_num = $.-1;
		$line=~/^(\d+)\s+/;
        @line = split/\s+/,$line;
        $chr = $line[0]; $pos = $line[1];
        push @chr,$chr;
        if ($hang_num == 1){
			$chr_star{$chr} = $hang_num;
		}
        $key = "$chr	$pos";
        $s4{$key} = $line[4]; $s5{$key} = $line[5]; $s6{$key} = $line[6]; $s7{$key} = $line[7]; 
    }
}
$chr_end{$chr} = $hang_num;
###########
 my (@array2,$i,$j,@array1);
open IN2,"$ARGV[1]" or die "can't open the input file 3";
while (<IN2>){
	chomp();
	if ($.>1){
		$j=$.-1;
		@array1=split/\s+/,$_; 
		foreach $i(0..$#array1){
			$array2[$j][$i]=$array1[$i];
		}
    }
}
$end{$chr}=$array2[$chr_end{$ARGV[3]}][1];
###########
my ($output,$sequence,$filename_in,@arrays,$window_num,@array,@source1);
foreach $sample_num($ARGV[4]..$ARGV[5]){
	print 'Staring sample'."$sample_num\n";
	$filename_in = 'step1_modify_'."$sample_num".'_'."$ARGV[3]".'.txt';
	$hang_in1=0;
	open IN4,"$filename_in" or die "can't open the file 4";
	@arrays = <IN4>;$window_num = $#arrays+1;my $hang_num = 0; 
	$window_num = $#arrays+1;my $hang_num = 0; 
	open IN4,"$filename_in" or die "can't open the file 4";
    while (<IN4>){ 
	    chomp();$hang_num++;
	   @array = split/\s+/,$_;
	   if ($hang_num==1){
			$win_num[$hang_num] = $array[1];
			$win_star[$hang_num] = 1;
			$win_end[$hang_num] = $array[1]*$ARGV[2];
			$source1[$hang_num] = $array[9];
			$source2[$hang_num] = $array[10];
		}elsif($hang_num>=2 and $hang_num < $window_num){
			$win_num[$hang_num] = $array[1];
            $win_star[$hang_num] = $win_num[$hang_num-1]*$ARGV[2]+1;
            $win_end[$hang_num] = $array[1]*$ARGV[2];
			$source1[$hang_num] = $array[9];
			$source2[$hang_num] = $array[10];
        }elsif($hang_num == $window_num){
    	    $win_num[$hang_num] = $array[1];
    	    $win_star[$hang_num] = $win_num[$hang_num-1]*$ARGV[2]+1;
    	    $win_end[$hang_num] = $end{$ARGV[3]};
			$source1[$hang_num] = $array[9];
			$source2[$hang_num] = $array[10];	
        } 
	}
	###########
    my @hang = ();my @snp = ();
    my ($hang,$mark_num_all,$mark_num_ture,$parent7,$parent6,$parent4,$parent5,$marker_use,$n);$parent4=0;$parent5=0;$parent6=0;$parent7=0;$mark_num_all=0;$mark_num_ture=0;
    $hang_num = 0;my $n = 1;
    my $filename_out = 'step2_modify_'."$sample_num".'_'."$ARGV[3]".'.txt';
	open OUT,">$filename_out";
	foreach $sequence($chr_star{$ARGV[3]}..$chr_end{$ARGV[3]}){
		$chr=$array2[$sequence][0];$hang_num++;
		$pos = $array2[$sequence][1];$key = "$chr	$pos";
		if ($pos >= $win_star[$n] and $pos <= $win_end[$n]){
			if ($source1[$n] eq "none" and $source2[$n] eq "none" ){
				$mark_num_all++;
				push @geno_hang,$sequence;
				if ($array2[$sequence][$sample_num+3] eq "H" or $array2[$sequence][$sample_num+3]=~/[A T G C]/){
					$mark_num_ture++;
				}
			}elsif ($source1[$n]=~/[A B]/ and $source2[$n]=~/[C D]/){
				$mark_num_all++;
				push @geno_hang,$sequence;
				if ($array2[$sequence][$sample_num+3] eq "H" or $array2[$sequence][$sample_num+3]=~/[A T G C]/){
					$mark_num_ture++;
				}
			}else{
				push @geno_hang,$sequence;
				$mark_num_all++;
				if ($array2[$sequence][$sample_num+3] eq "H"){
					$mark_num_ture++;
					if ($source1[$n] eq "A" and $s6{$key} ne $s7{$key} and $s4{$key} ne $s7{$key}){
						$parent7++;
					}elsif($source1[$n] eq "A" and $s6{$key} ne $s7{$key} and $s4{$key} ne $s6{$key}){
						$parent6++;
					}elsif($source1[$n] eq "B" and $s6{$key} ne $s7{$key} and $s5{$key} ne $s7{$key}){
						$parent7++;
					}elsif($source1[$n] eq "B" and $s6{$key} ne $s7{$key} and $s5{$key} ne $s6{$key}){
						$parent6++;
					}elsif($source2[$n] eq "C" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s4{$key}){
						$parent4++;
					}elsif($source2[$n] eq "C" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s5{$key}){
						$parent5++;
					}elsif($source2[$n] eq "D" and $s4{$key} ne $s5{$key} and $s7{$key} ne $s4{$key}){
						$parent4++;
					}elsif($source2[$n] eq "D" and $s4{$key} ne $s5{$key} and $s7{$key} ne $s5{$key}){
						$parent5++;
					}
				}elsif($array2[$sequence][$sample_num+3]=~/[A T G C]/){
					$mark_num_ture++;
					if($source1[$n] eq "A" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s7{$key} and $array2[$sequence][$sample_num+3] eq $s6{$key}){
						$parent6++;
					}elsif($source1[$n] eq "A" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s6{$key} and $array2[$sequence][$sample_num+3] eq $s7{$key}){
						$parent7++;
					}elsif($source1[$n] eq "B" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s7{$key} and $array2[$sequence][$sample_num+3] eq $s7{$key}){
						$parent7++;
					}elsif($source1[$n] eq "B" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s6{$key} and $array2[$sequence][$sample_num+3] eq $s6{$key}){
						$parent6++;
					}elsif($source2[$n] eq "C" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s6{$key} and $array2[$sequence][$sample_num+3] eq $s5{$key}){
						$parent5++;
					}elsif($source2[$n] eq "C" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s7{$key} and $array2[$sequence][$sample_num+3] eq $s4{$key}){
						$parent4++;
					}elsif($source2[$n] eq "D" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s7{$key} and $array2[$sequence][$sample_num+3] eq $s5{$key}){
						$parent5++;
					}elsif($source2[$n] eq "D" and $s4{$key} ne $s5{$key} and $s6{$key} ne $s7{$key} and $s4{$key} eq $s6{$key} and $array2[$sequence][$sample_num+3] eq $s4{$key}){
						$parent4++;
					}elsif($s4{$key} eq $s5{$key} and $s4{$key} eq $s6{$key} and $s4{$key} ne $s7{$key} and $array2[$sequence][$sample_num+3] eq $s7{$key} and $source2[$n] ne "C" and $source2[$n] ne "D"){
						$parent7++;
					}elsif($s4{$key} eq $s5{$key} and $s4{$key} eq $s7{$key} and $s4{$key} ne $s6{$key} and $array2[$sequence][$sample_num+3] eq $s6{$key} and $source2[$n] ne "C" and $source2[$n] ne "D"){
						$parent6++;
					}elsif($s4{$key} eq $s6{$key} and $s4{$key} eq $s7{$key} and $s4{$key} ne $s5{$key} and $array2[$sequence][$sample_num+3] eq $s5{$key} and $source1[$n] ne "A" and $source1[$n] ne "B"){
						$parent5++;
					}elsif($s4{$key} ne $s5{$key} and $s5{$key} eq $s7{$key} and $s5{$key} eq $s7{$key} and $array2[$sequence][$sample_num+3] eq $s4{$key} and $source1[$n] ne "A" and $source1[$n] ne "B"){
						$parent4++;
					}
				}
			}
			if ($sequence == $chr_end{$ARGV[3]}){
				$marker_use = $parent4+$parent5+$parent6+$parent7;
				$snp[$n][0] = $ARGV[3];$snp[$n][1]=$win_star[$n];$snp[$n][2]=$win_end[$n];$snp[$n][3]=$mark_num_all;$snp[$n][4]=$mark_num_ture;
				$snp[$n][5]=$marker_use;$snp[$n][6]=$parent4;$snp[$n][7]=$parent5;$snp[$n][8]=$parent6;$snp[$n][9]=$parent7;
				$snp[$n][12]=$geno_hang[0];$snp[$n][13]=$geno_hang[$#geno_hang];
				if ($source1[$n] eq "none" and $source2[$n] eq "none" ){
					$snp[$n][10]=$source1[$n];$snp[$n][11]=$source2[$n];
				}elsif ($source1[$n]=~/[A B]/ and $source2[$n]=~/[C D]/){
					$snp[$n][10]=$source1[$n];$snp[$n][11]=$source2[$n];
				}else{
					if($source1[$n] eq "none"){
						$snp[$n][11]=$source2[$n];
						if(max($parent4,$parent5)==$parent4 and $parent4 >= 15){
							$snp[$n][10]="A";
						}elsif(max($parent4,$parent5)==$parent5 and $parent5 >= 15){
							$snp[$n][10]="B";
						}else{
							$snp[$n][10]="none";
						}
					}elsif($source2[$n] eq "none"){
						$snp[$n][10]=$source1[$n];
						if(max($parent6,$parent7)==$parent6 and $parent6 >= 15){
							$snp[$n][11]="C";
						}elsif(max($parent6,$parent7)==$parent7 and $parent7 >= 15){
							$snp[$n][11]="D";
						}else{
							$snp[$n][11]="none";
						}
					}
				}	
			}
		}elsif($pos>$win_end[$n]){
	 		$marker_use = $parent4+$parent5+$parent6+$parent7;
			$snp[$n][0] = $ARGV[3];$snp[$n][1]=$win_star[$n];$snp[$n][2]=$win_end[$n];$snp[$n][3]=$mark_num_all;$snp[$n][4]=$mark_num_ture;
			$snp[$n][5]=$marker_use;$snp[$n][6]=$parent4;$snp[$n][7]=$parent5;$snp[$n][8]=$parent6;$snp[$n][9]=$parent7;$snp[$n][12]=$geno_hang[0];$snp[$n][13]=$geno_hang[$#geno_hang];
			if ($source1[$n] eq "none" and $source2[$n] eq "none" ){
				$snp[$n][10]=$source1[$n];$snp[$n][11]=$source2[$n];
			}elsif ($source1[$n]=~/[A B]/ and $source2[$n]=~/[C D]/){
				$snp[$n][10]=$source1[$n];$snp[$n][11]=$source2[$n];
			}else{
				if($source1[$n] eq "none"){
					$snp[$n][11]=$source2[$n];
					if(max($parent4,$parent5)==$parent4 and $parent4 >= 15){
						$snp[$n][10]="A";
					}elsif(max($parent4,$parent5)==$parent5 and $parent5 >= 15){
						$snp[$n][10]="B";
					}else{
						$snp[$n][10]="none";
					}
				}elsif($source2[$n] eq "none"){
					$snp[$n][10]=$source1[$n];
					if(max($parent6,$parent7)==$parent6 and $parent6 >= 15){
						$snp[$n][11]="C";
					}elsif(max($parent6,$parent7)==$parent7 and $parent7 >= 15){
						$snp[$n][11]="D";
					}else{
						$snp[$n][11]="none";
					}
				}
			}	
	 		$n++;$parent4=0;$parent5=0;$parent6=0;$parent7=0;$mark_num_all=0;$mark_num_ture=0;@geno_hang=();
	 		redo;
        }
	}
	#填充头缺失
	$num1=0;
	foreach $window_num(1..$n){
		if ($snp[1][10] eq "none"){
			if ($snp[$window_num][10] eq "none"){
				$num1++;
			}else{
				foreach (1..$num1){
					$snp[$_][10]=$snp[$num1+1][10];
				}
				last;
			}
		}
	}
	$num1=0;
	foreach $window_num(1..$n){
		if ($snp[1][11] eq "none"){
			if ($snp[$window_num][11] eq "none"){
				$num1++;
			}else{
				foreach (1..$num1){
					$snp[$_][11]=$snp[$num1+1][11];
				}
				last;
			}
		}
	}
	@n=sort { $b <=> $a } (1..$n);
	$num1=0;
	foreach $window_num(@n){
		if ($snp[$n[0]][10] eq "none"){
			if ($snp[$window_num][10] eq "none"){
				$num1++;
			}else{
				foreach ($n[0]-$num1+1..$n[0]){
					$snp[$_][10]=$snp[$n[0]-$num1][10];
				}
				last;
			}
		}
	}
	$num1=0;
	foreach $window_num(@n){
		if ($snp[$n[0]][11] eq "none"){
			if ($snp[$window_num][11] eq "none"){
				$num1++;
			}else{
				foreach ($n[0]-$num1+1..$n[0]){
					$snp[$_][11]=$snp[$n[0]-$num1][11];
				}
				last;
			}
		}
	}
	$num1=0;
	foreach $window_num(1..$n){
		if ($snp[$window_num][10] eq "none"){
			$num1++;
		}else{
			if ($window_num-$num1==1){
			}elsif($snp[$window_num-1][10] eq "none"){
				if ($snp[$window_num][10] eq $snp[$window_num-$num1-1][10]){
					foreach ($window_num-$num1..$window_num-1){
						$snp[$_][10]="$snp[$window_num][10]";
					}
				}else{
					if (($num1%2)==0){
						foreach ($window_num-$num1..$window_num-$num1+int($num1/2)-1){
							$snp[$_][10]="$snp[$window_num-$num1-1][10]";
						}
						foreach ($window_num-$num1+int($num1/2)..$window_num-1){
							$snp[$_][10]="$snp[$window_num][10]";
						}
					}else{
						foreach ($window_num-$num1..$window_num-$num1+int($num1/2)){
							$snp[$_][10]="$snp[$window_num-$num1-1][10]";
						}
						foreach ($window_num-$num1+int($num1/2)+1..$window_num-1){
							$snp[$_][10]= "$snp[$window_num][10]";
						}
					}
				}
			}
			$num1=0;
		}
	}
	$num1=0;
	foreach $window_num(1..$n){
		if ($snp[$window_num][11] eq "none"){
			$num1++;
		}else{
			if ($window_num-$num1==1){
			}elsif($snp[$window_num-1][11] eq "none"){
				if ($snp[$window_num][11] eq $snp[$window_num-$num1-1][11]){
					foreach ($window_num-$num1..$window_num-1){
						$snp[$_][11]="$snp[$window_num][11]";
					}
				}else{
					if (($num1%2)==0){
						foreach ($window_num-$num1..$window_num-$num1+int($num1/2)-1){
							$snp[$_][11]="$snp[$window_num-$num1-1][11]";
						}
						foreach ($window_num-$num1+int($num1/2)..$window_num-1){
							$snp[$_][11]="$snp[$window_num][11]";
						}
					}else{
						foreach ($window_num-$num1..$window_num-$num1+int($num1/2)){
							$snp[$_][11]="$snp[$window_num-$num1-1][11]";
						}
						foreach ($window_num-$num1+int($num1/2)+1..$window_num-1){
							$snp[$_][11]= "$snp[$window_num][11]";
						}
					}
				}
			}
			$num1=0;
		}
	}
	
    ######
	my (@parent);
	foreach $window_num(1..$n){print OUT "$snp[$window_num][0]	$snp[$window_num][1]	$snp[$window_num][2]	$snp[$window_num][3]	$snp[$window_num][4]	$snp[$window_num][5]	$snp[$window_num][6]	$snp[$window_num][7]	$snp[$window_num][8]	$snp[$window_num][9]	$snp[$window_num][10]	$snp[$window_num][11]	$snp[$window_num][12]	$snp[$window_num][13]\n";
	}	
	########
	my ($staring,$ending,$origin1,$origin2,$position,$file);
	$file = "$ARGV[3]".'_'."$sample_num".'modify';
	open OUT1,">$file.parent1";
	open OUT2,">$file.parent2";
	foreach $window_num(1..$n){
		if($snp[$window_num][0] == $ARGV[3]){
			if($window_num<$n){
				$staring = int($snp[$window_num][1]/$ARGV[2]);$ending = int($snp[$window_num][2]/$ARGV[2]);
			}elsif($window_num == $n){
				$staring = int($snp[$window_num][1]/$ARGV[2]);$ending = int($snp[$window_num][2]/$ARGV[2])+1;
				}
			if ($ending - $staring == 1){
				$position = $ending;
				print OUT1 "$position\t$snp[$window_num][10]\n";
				print OUT2 "$position\t$snp[$window_num][11]\n";
			}else{
				$origin1 = $snp[$window_num][10];$origin2 = $snp[$window_num][11];
				foreach ($staring+1..$ending){
					$position = $_ ;
					print OUT1 "$position\t$origin1\n";
					print OUT2 "$position\t$origin2\n";
				}
			}
		}
	}
	close OUT1;
	close OUT2;     
}
#######
