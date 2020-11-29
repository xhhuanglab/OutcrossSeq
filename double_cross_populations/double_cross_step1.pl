my ($line,@line,$hang_num,$chr,@chr,$key,$pos);our (%s4,%s5,%s6,%s7,%chr_star,%chr_end);
open IN1,"$ARGV[0]" or die "can't open the input file 1";
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
close IN1;
###########
my %all_windows_num;my %seq;
my ($window_length,$all_windows_num,$array,$mark_num_all,$mark_num_ture);
my $windows_num=1;my $s4_num=0;my $s5_num=0;my $s6_num=0;my $s7_num=0;my $marker_use=0;
$window_length = $ARGV[2];
if ($pos%$window_length == 0){
	$all_windows_num{$ARGV[3]} = int ($pos/$window_length);   
}else{
    $all_windows_num{$ARGV[3]} = int ($pos/$window_length)+1;
}
#######
 my (@array,$i,$j,@array1);
open IN2,"$ARGV[1]" or die "can't open the input file 3";
while (<IN2>){
	chomp();
	if ($.>1){
		$j=$.-1;
		@array1=split/\s+/,$_; 
		foreach $i(0..$#array1){
			$array[$j][$i]=$array1[$i];
		}
    }
}
if (int($array[$chr_end{$ARGV[3]}][1]/$window_length)+1 != $all_windows_num{$ARGV[3]}){
	$all_windows_num{$ARGV[3]}=int($array[$chr_end{$ARGV[3]}][1]/$window_length)+1;
}
###########
###########
print "Begaining\n";
my ($output,$sequence);
foreach $sample_num($ARGV[4]..$ARGV[5]){
	print 'Staring sample'."$sample_num\n";
	$output = 'step1_modify_'."$sample_num".'_'."$ARGV[3]";
	open OUT,">$output.txt";
	my @hang = ();my @markers = (); $hang[0][0]="Windows_num";
	$hang[0][1]="Mark_num_all";$hang[0][2]="Mark_num_ture";$hang[0][3]="marker_use";
	$hang[0][4]="S4";$hang[0][5]="S5";$hang[0][6]="S6";$hang[0][7]="S7";
	my $windows_num=1;my $s4_num=0;my $s5_num=0;my $s6_num=0;my $s7_num=0;my $marker_use=0;my $mark_num_all=0;my $mark_num_ture=0;
	foreach $sequence($chr_star{$ARGV[3]}..$chr_end{$ARGV[3]}){    
	    $pos = $array[$sequence][1];$chr=$array[$sequence][0];
	    $key = "$chr	$pos";
	    if (int($pos/$window_length)+1 == $windows_num and $windows_num <= $all_windows_num{$ARGV[3]}){
			$mark_num_all++;
			if ($array[$sequence][$sample_num+3] eq "H"){
				$mark_num_ture++;
	            if ($s4{$key} eq $s5{$key} and $s4{$key} eq $s6{$key} and $s4{$key} ne $s7{$key}){
					$s7_num++;
                }elsif ($s4{$key} eq $s5{$key} and $s4{$key} eq $s7{$key} and $s4{$key} ne $s6{$key}){
					$s6_num++; 
                }elsif ($s4{$key} eq $s6{$key} and $s4{$key} eq $s7{$key} and $s4{$key} ne $s5{$key}){
					$s5_num++; 
                }elsif ($s5{$key} eq $s6{$key} and $s5{$key} eq $s7{$key} and $s4{$key} ne $s5{$key}){
					$s4_num++; 
                }
            }elsif ($array[$sequence][$sample_num+3]=~/[A T G C]/){
            	$mark_num_ture++;
            	if ($s4{$key} eq $s5{$key} and $s4{$key} eq $s6{$key} and $s4{$key} ne $s7{$key}){
            		if ($array[$sequence][$sample_num+3] eq $s7{$key}){
						$s7_num++;
					}
            	 }elsif ($s4{$key} eq $s5{$key} and $s4{$key} eq $s7{$key} and $s4{$key} ne $s6{$key}){
            	    if ($array[$sequence][$sample_num+3] eq $s6{$key}){	
						$s6_num++;
					}
            	 }elsif ($s4{$key} eq $s6{$key} and $s4{$key} eq $s7{$key} and $s4{$key} ne $s5{$key}){
            	 	if ($array[$sequence][$sample_num+3] eq $s5{$key}){
						$s5_num++;
					}
                 }elsif ($s5{$key} eq $s6{$key} and $s5{$key} eq $s7{$key} and $s4{$key} ne $s5{$key}){
            		if ($array[$sequence][$sample_num+3] eq $s4{$key}){
						$s4_num++;
					}
                 }
            }
            if($sequence == $chr_end{$ARGV[3]}  and int($array[$sequence][1]/$window_length)+1 == $all_windows_num{$ARGV[3]}){
               $marker_use = $s4_num+$s5_num+$s6_num+$s7_num;$hang[$all_windows_num{$ARGV[3]}][0]=$all_windows_num{$ARGV[3]};
			   $hang[$all_windows_num{$ARGV[3]}][1]=$mark_num_all;$hang[$all_windows_num{$ARGV[3]}][2]=$mark_num_ture;
			   $hang[$all_windows_num{$ARGV[3]}][3]=$marker_use;$hang[$all_windows_num{$ARGV[3]}][4]=$s4_num;
			   $hang[$all_windows_num{$ARGV[3]}][5]=$s5_num;$hang[$all_windows_num{$ARGV[3]}][6]=$s6_num;
			   $hang[$all_windows_num{$ARGV[3]}][7]=$s7_num;
           }
        }elsif(int($array[$sequence][1]/$window_length)+1 != $windows_num and int($array[$sequence][1]/$window_length)+1 <= $all_windows_num{$ARGV[3]}){
               $marker_use = $s4_num+$s5_num+$s6_num+$s7_num;$hang[$windows_num][0]=$windows_num;
			   $hang[$windows_num][1]=$mark_num_all;$hang[$windows_num][2]=$mark_num_ture;
			   $hang[$windows_num][3]=$marker_use;$hang[$windows_num][4]=$s4_num;
			   $hang[$windows_num][5]=$s5_num;$hang[$windows_num][6]=$s6_num;$hang[$windows_num][7]=$s7_num; 
		       $windows_num++;$mark_num_all = 0;$mark_num_ture = 0;$s4_num=0;$s5_num=0;$s6_num=0;$s7_num=0;
			  
               redo;
        }
	}
	my ($mark_windows,$n,$nomiss,$allsite,$snp_marker,$S4_num,$S5_num,$S6_num,$S7_num,@markers); 
	$markers[0][8]="Windows_num";$markers[0][0]="End_window_num";$markers[0][1]="Mark_num_all";$markers[0][2]="Mark_num_ture";$markers[0][3]="Marker_use";$markers[0][4]="S4";$markers[0][5]="S5";$markers[0][6]="S6";$markers[0][7]="S7";	$markers[0][9]="p0";
	foreach $mark_windows(1..$windows_num){
		if ($mark_windows==1 and $hang[$mark_windows][2]>=150){
			$n=1;$markers[$n][8]=$n;
			foreach(0..7){
			   $markers[$n][$_]=$hang[$mark_windows][$_];
			   }
		}elsif ($mark_windows > 1 and $mark_windows <  $windows_num and $hang[$mark_windows][2] >= 150 and $hang[$mark_windows-1][2] >= 150 ){
			$n++;$markers[$n][8]=$n;
			foreach(0..7){
			   $markers[$n][$_]=$hang[$mark_windows][$_];
			}
		}else{ 
			$nomiss+=$hang[$mark_windows][1];
			$allsite+=$hang[$mark_windows][2];
			$snp_marker+=$hang[$mark_windows][3];
			$S4_num+=$hang[$mark_windows][4];
			$S5_num+=$hang[$mark_windows][5];
			$S6_num+=$hang[$mark_windows][6];
			$S7_num+=$hang[$mark_windows][7];
			if ($allsite >= 150 and $mark_windows < $windows_num){
				 $n++;$markers[$n][8]=$n;
				 $markers[$n][0]=$hang[$mark_windows][0];$markers[$n][1]=$nomiss;$markers[$n][2]=$allsite;$markers[$n][3]=$snp_marker;$markers[$n][4]=$S4_num;$markers[$n][5]=$S5_num;$markers[$n][6]=$S6_num;$markers[$n][7]=$S7_num;                  
				 $S7_num=0;$S6_num=0;$S5_num=0;$S4_num=0;$allsite=0;$nomiss=0;$snp_marker=0;
			}
			if($mark_windows==$windows_num){
				$n++;$markers[$n][8]=$n;
				$markers[$n][0]=$hang[$mark_windows][0];$markers[$n][1]=$nomiss;$markers[$n][2]=$allsite;$markers[$n][3]=$snp_marker;$markers[$n][4]=$S4_num;$markers[$n][5]=$S5_num;$markers[$n][6]=$S6_num;$markers[$n][7]=$S7_num;
			}
		}
		
	}
	my ($source,@parent,@end_star);
	foreach (1..$n){
		$source1 = "none";$source2 = "none";
		if ($markers[$_][2] > 0 and $markers[$_][3]/$markers[$_][2]>=0.08){
			if (($markers[$_][4]+$markers[$_][5])/$markers[$_][3] >= 0.4 or ($markers[$_][4]+$markers[$_][5])>=20){
				if($markers[$_][4]/($markers[$_][4]+$markers[$_][5])>=0.8){
					$source1="A";
				}elsif($markers[$_][5]/($markers[$_][4]+$markers[$_][5])>=0.8){
					$source1="B";
				}elsif($markers[$_][5]/($markers[$_][4]+$markers[$_][5])>=0.4 and $markers[$_][5]/($markers[$_][4]+$markers[$_][5])<=0.6){
					$source1="none";
				}
			}
			if (($markers[$_][6]+$markers[$_][7])/$markers[$_][3] >= 0.4 or ($markers[$_][6]+$markers[$_][7])>=20){
				if($markers[$_][6]/($markers[$_][6]+$markers[$_][7])>=0.8){
					$source2="C";
				}elsif($markers[$_][7]/($markers[$_][6]+$markers[$_][7])>=0.8){
					$source2="D";
				}elsif($markers[$_][6]/($markers[$_][6]+$markers[$_][7])>=0.4 and $markers[$_][6]/($markers[$_][6]+$markers[$_][7])<=0.6){
					$source2="none";
				}
			}
		}
				
		$markers[$_][9]="$source1";
		$markers[$_][10]="$source2"; 
	}
	@star=();$star_num=0;$star_num_list=0;
	foreach $i(1..$n){
		if ($i==1){
			$star=$markers[$i][9];
		}
		if ($markers[$i][9] eq $star){
			$star_num++;
		}else{
			$star_num_list++;	$star=$markers[$i][9];
			if($i-$star_num==1){
				$star[$star_num_list][0]=0;
			}else{
				$star[$star_num_list][0]=$markers[$i-$star_num-1][0];
			}
			$star[$star_num_list][1]=$markers[$i-1][0];
			$star[$star_num_list][2]=$markers[$i-1][0]-$markers[$i-$star_num-1][0];
			$star[$star_num_list][3]=$markers[$i-1][9];
			$star[$star_num_list][4]=$markers[$i-$star_num][8];
			$star[$star_num_list][5]=$markers[$i-1][8];
			$star_num=1;
		}
		if ($i==$n){
			$star_num_list++;				
			$star[$star_num_list][0]=$markers[$i-$star_num][0];
			$star[$star_num_list][1]=$markers[$i][0];
			$star[$star_num_list][2]=$markers[$i][0]-$markers[$i-$star_num][0];
			$star[$star_num_list][3]=$markers[$i][9];
			$star[$star_num_list][4]=$markers[$i-$star_num][8];
			$star[$star_num_list][5]=$i;
		}
	}	
	@star_num_list_geno=();
	@star_num_list_none=();
	foreach $i(1..$star_num_list){
		if ($star[$i][3]=~/[A B]/){
			push @star_num_list_geno,$i;
		}else{
			push @star_num_list_none,$i;
		}
	}
	foreach $i(0..$#star_num_list_geno){
		if ($star[$star_num_list_geno[$i]][0]==0){
			if ($star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i+1]][3]){
				if ($star[$star_num_list_geno[$i]][2]*$window_length<=300000 or $star[$star_num_list_geno[$i]][4]==$star[$star_num_list_geno[$i]][5]){
					$num1=0;
					foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
						if ($markers[$_][4]<20 and $markers[$_][5]<20){
							$num1++;
						}
					}
					if ($num1==($star[$star_num_list_geno[$i]][5]-$star[$star_num_list_geno[$i]][4]+1)){
						foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
							$markers[$_][9]="none";
						}
					}
				}
			}
		}elsif($star[$star_num_list_geno[$i]][5]==$n){
			if ($star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i-1]][3]){
				if($star[$star_num_list_geno[$i]][2]*$window_length<=300000 or $star[$star_num_list_geno[$i]][4]==$star[$star_num_list_geno[$i]][5]){
					$num1=0;
					foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
						if ($markers[$_][4]<20 and $markers[$_][5]<20){
							$num1++;
						}
					}
					if ($num1==($star[$star_num_list_geno[$i]][5]-$star[$star_num_list_geno[$i]][4]+1)){
						foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
							$markers[$_][9]="none";
						}
					}
				}
			}
		}else{
			if($star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i-1]][3] and $star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i+1]][3]){
				if($star[$star_num_list_geno[$i]][2]*$window_length<=300000 or $star[$star_num_list_geno[$i]][4]==$star[$star_num_list_geno[$i]][5]){
					$num1=0;
					foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
						if ($markers[$_][4]<20 and $markers[$_][5]<20){
							$num1++;
						}
					}
					if ($num1==($star[$star_num_list_geno[$i]][5]-$star[$star_num_list_geno[$i]][4]+1)){
						foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
							$markers[$_][9]="none";
						}
					}
				}
			}
		}
	}
	foreach $i(0..$#star_num_list_none){
		if ($star[$star_num_list_none[$i]][0]==0){
			if($star[$star_num_list_none[$i]][2]*$window_length<=300000 or $star[$star_num_list_none[$i]][4]==$star[$star_num_list_none[$i]][5]){
				foreach ($star[$star_num_list_none[$i]][4]..$star[$star_num_list_none[$i]][5]){
					$markers[$_][9]=$star[$star_num_list_geno[0]][3];
				}
			}
		}elsif($star[$star_num_list_none[$i]][5]==$n){
			if($star[$star_num_list_none[$i]][2]*$window_length<=300000 or $star[$star_num_list_none[$i]][4]==$star[$star_num_list_none[$i]][5]){
				foreach ($star[$star_num_list_none[$i]][4]..$star[$star_num_list_none[$i]][5]){
					$markers[$_][9]=$star[$star_num_list_geno[$#star_num_list_geno]][3];
				}
			}
		}else{
			if ($markers[$star[$star_num_list_none[$i]][4]-1][9] eq $markers[$star[$star_num_list_none[$i]][5]+1][9]){
				if($star[$star_num_list_none[$i]][2]*$window_length<=300000 or $star[$star_num_list_none[$i]][4]==$star[$star_num_list_none[$i]][5]){
					foreach ($star[$star_num_list_none[$i]][4]..$star[$star_num_list_none[$i]][5]){
						$markers[$_][9]=$markers[$star[$star_num_list_none[$i]][5]+1][9];
					}
				}
			}
		}
	}
	@star=();$star_num=0;$star_num_list=0;
	foreach $i(1..$n){
		if ($i==1){
			$star=$markers[$i][10];
		}
		if ($markers[$i][10] eq $star){
			$star_num++;
		}else{
			$star_num_list++;	$star=$markers[$i][10];
			if($i-$star_num==1){
				$star[$star_num_list][0]=0;
			}else{
				$star[$star_num_list][0]=$markers[$i-$star_num-1][0];
			}
			$star[$star_num_list][1]=$markers[$i-1][0];
			$star[$star_num_list][2]=$markers[$i-1][0]-$markers[$i-$star_num-1][0];
			$star[$star_num_list][3]=$markers[$i-1][10];
			$star[$star_num_list][4]=$markers[$i-$star_num][8];
			$star[$star_num_list][5]=$markers[$i-1][8];
			$star_num=1;
		}
		if ($i==$n){
			$star_num_list++;				
			$star[$star_num_list][0]=$markers[$i-$star_num][0];
			$star[$star_num_list][1]=$markers[$i][0];
			$star[$star_num_list][2]=$markers[$i][0]-$markers[$i-$star_num][0];
			$star[$star_num_list][3]=$markers[$i][10];
			$star[$star_num_list][4]=$markers[$i-$star_num][8];
			$star[$star_num_list][5]=$i;
		}
	}	
	@star_num_list_geno=();
	@star_num_list_none=();
	foreach $i(1..$star_num_list){	
		if ($star[$i][3]=~/[C D]/){
			push @star_num_list_geno,$i;
		}else{
			push @star_num_list_none,$i;
		}
	}
	foreach $i(0..$#star_num_list_geno){
		if ($star[$star_num_list_geno[$i]][0]==0){
			if ($star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i+1]][3]){
				if ($star[$star_num_list_geno[$i]][2]*$window_length<=300000 or $star[$star_num_list_geno[$i]][4]==$star[$star_num_list_geno[$i]][5]){
					$num1=0;
					foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
						if ($markers[$_][6]<20 and $markers[$_][7]<20){
							$num1++;
						}
					}
					if ($num1==($star[$star_num_list_geno[$i]][5]-$star[$star_num_list_geno[$i]][4]+1)){
						foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
							$markers[$_][10]="none";
						}
					}
				}
			}
		}elsif($star[$star_num_list_geno[$i]][5]==$n){
			if ($star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i-1]][3]){
				if($star[$star_num_list_geno[$i]][2]*$window_length<=300000 or $star[$star_num_list_geno[$i]][4]==$star[$star_num_list_geno[$i]][5]){
					$num1=0;
					foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
						if ($markers[$_][6]<20 and $markers[$_][7]<20){
							$num1++;
						}
					}
					if ($num1==($star[$star_num_list_geno[$i]][5]-$star[$star_num_list_geno[$i]][4]+1)){
						foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
							$markers[$_][10]="none";
						}
					}
				}
			}
		}else{
			if($star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i-1]][3] and $star[$star_num_list_geno[$i]][3] ne $star[$star_num_list_geno[$i+1]][3]){
				if($star[$star_num_list_geno[$i]][2]*$window_length<=300000 or $star[$star_num_list_geno[$i]][4]==$star[$star_num_list_geno[$i]][5]){
					$num1=0;
					foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
						if ($markers[$_][6]<20 and $markers[$_][7]<20){
							$num1++;
						}
					}
					if ($num1==($star[$star_num_list_geno[$i]][5]-$star[$star_num_list_geno[$i]][4]+1)){
						foreach ($star[$star_num_list_geno[$i]][4]..$star[$star_num_list_geno[$i]][5]){
							$markers[$_][10]="none";
						}
					}
				}
			}
		}
	}
	
	foreach $i(0..$#star_num_list_none){
		if ($star[$star_num_list_none[$i]][0]==0){
			if($star[$star_num_list_none[$i]][2]*$window_length<=300000 or $star[$star_num_list_none[$i]][4]==$star[$star_num_list_none[$i]][5]){
				foreach ($star[$star_num_list_none[$i]][4]..$star[$star_num_list_none[$i]][5]){
					$markers[$_][10]=$star[$star_num_list_geno[0]][3];
				}
			}
		}elsif($star[$star_num_list_none[$i]][5]==$n){
			if($star[$star_num_list_none[$i]][2]*$window_length<=300000 or $star[$star_num_list_none[$i]][4]==$star[$star_num_list_none[$i]][5]){
				foreach ($star[$star_num_list_none[$i]][4]..$star[$star_num_list_none[$i]][5]){
					$markers[$_][10]=$star[$star_num_list_geno[$#star_num_list_geno]][3];
				}
			}
		}else{
			if ($markers[$star[$star_num_list_none[$i]][4]-1][10] eq $markers[$star[$star_num_list_none[$i]][5]+1][10]){
				if($star[$star_num_list_none[$i]][2]*$window_length<=300000 or $star[$star_num_list_none[$i]][4]==$star[$star_num_list_none[$i]][5]){
					foreach ($star[$star_num_list_none[$i]][4]..$star[$star_num_list_none[$i]][5]){
						$markers[$_][10]=$markers[$star[$star_num_list_none[$i]][4]-1][10];
					}
				}
			}
		}
	}
	
	$num1=0;$num2=0;
	foreach $i(1..$n){
		if ($markers[1][9] eq "none" and $markers[1][10] eq "none"){
			if ($markers[$i][9] eq "none"){$num1++;}else{last;}
			
		}
	}
	foreach $i(1..$n){
		if ($markers[1][9] eq "none" and $markers[1][10] eq "none"){
			if ($markers[$i][10] eq "none"){$num2++;}else{last;}
		}
	}
	if ($num1>=$num2){
		foreach (1..$num2){
			$markers[$_][10]=$markers[$num2+1][10];
		}
	}else{
		foreach (1..$num1){
			$markers[$_][9]=$markers[$num1+1][9];
		}
	}
	@n=sort { $b <=> $a } (1..$n);
	$num1=0;$num2=0;
	foreach $i(@n){
		if ($markers[$n[0]][9] eq "none" and $markers[$n[0]][10] eq "none"){
			if ($markers[$i][9] eq "none"){$num1++;}else{last;}
		}
	}
	foreach $i(@n){
		if ($markers[$n[0]][9] eq "none" and $markers[$n[0]][10] eq "none"){
			if ($markers[$i][10] eq "none"){$num2++;}else{last;}
		}
	}
	if ($num1>=$num2){
		foreach ($n[0]-$num2+1..$n[0]){
			$markers[$_][10]=$markers[$n[0]-$num2][10];
		}
	}else{
		foreach ($n[0]-$num1+1..$n[0]){
			$markers[$_][9]=$markers[$n[0]-$num1][9];
		}
	}	
	$num1=0;$num2_1=0;
	foreach $i(1..$n){
		if ($markers[$i][9] eq "none"){
			$num1++;
		}else{
			if ($i-$num1==1){
			}elsif($markers[$i-1][9] eq "none"){
				$num2_1++;
				$aa[$num2_1][0]=$markers[$i-$num1][8];$aa[$num2_1][1]=$markers[$i-1][8];
				$aa[$num2_1][2]=$markers[$i-$num1-1][9];$aa[$num2_1][3]=$markers[$i][9];
				$num1=0;
			}
			$num1=0;
		}
	}
	$num1=0;$num2_2=0;
	foreach $i(1..$n){
		if ($markers[$i][10] eq "none"){
			$num1++;
		}else{
			if ($i-$num1==1){
			}elsif($markers[$i-1][10] eq "none"){
				$num2_2++;
				$bb[$num2_2][0]=$markers[$i-$num1][8];$bb[$num2_2][1]=$markers[$i-1][8];
				$bb[$num2_2][2]=$markers[$i-$num1-1][10];$bb[$num2_2][3]=$markers[$i][10];
				$num1=0;
			}
			$num1=0;
		}
	}
	$num1=0;$num2_3=0;
	foreach $i(1..$n){
		if ($markers[$i][9] eq "none" and $markers[$i][10] eq "none"){
			$num1++;
		}else{
			if ($i-$num1==1){
			}elsif($markers[$i-1][9] eq "none" and $markers[$i-1][10] eq "none"){
				$num2_3++;
				$cc[$num2_3][0]=$markers[$i-$num1][8];$cc[$num2_3][1]=$markers[$i-1][8];
				$cc[$num2_3][2]=$markers[$i-$num1-1][9];$cc[$num2_3][3]=$markers[$i][9];
				$cc[$num2_3][4]=$markers[$i-$num1-1][10];$cc[$num2_3][5]=$markers[$i][10];
			}
			$num1=0;
		}
	}
	$num2=0;$num1=0;
	foreach $i(1..$num2_3){
		foreach $j(1..$num2_1){
			if ($cc[$i][0]>=$aa[$j][0] and $cc[$i][1]<=$aa[$j][1]){
				if ($aa[$j][2] eq $aa[$j][3]){
					$num1=$aa[$j][1]-$aa[$j][0];
					$geno1=$aa[$j][2];
				}else{
					$num1=0;
				}
			}
		}
		foreach $j(1..$num2_2){
			if ($cc[$i][0]>=$bb[$j][0] and $cc[$i][1]<=$bb[$j][1]){
				if ($bb[$j][2] eq $bb[$j][3]){
					$num2=$bb[$j][1]-$bb[$j][0];
					$geno2=$bb[$j][3];
				}else{
					$num2=0;
				}
			}
		}
		if ($num1>$num2 and $num2>0){
			foreach ($cc[$i][0]..$cc[$i][1]){
				$markers[$_][10]=$geno2;
			}
		}elsif($num2>$num1 and $num1>0){
			foreach ($cc[$i][0]..$cc[$i][1]){
				$markers[$_][9]=$geno1;
			}
		}elsif( $num1>0 and $num2==0){
			foreach ($cc[$i][0]..$cc[$i][1]){
				$markers[$_][9]=$geno1;
			}
		}elsif( $num2>0 and $num1==0){
			foreach ($cc[$i][0]..$cc[$i][1]){
				$markers[$_][10]=$geno2;
			}
		}
	}
	foreach $i(1..$n){print OUT "$markers[$i][8]	$markers[$i][0]	$markers[$i][1]	$markers[$i][2]	$markers[$i][3]	$markers[$i][4]	$markers[$i][5]	$markers[$i][6]	$markers[$i][7]	$markers[$i][9]	$markers[$i][10]\n";}
}
close OUT;
