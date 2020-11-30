#!/usr/bin/perl-w
$chr = $ARGV[0];
use List::Util qw/max min/;
open IN,"sample.list";
chomp(@array=<IN>); 
$sample_num = $#array+1;
$min_sample_num =int( $sample_num/4*0.7);
$max_sample_num =int( $sample_num/4*1.3);
open IN1,"Chr$chr.genotype" or die "can't open the input file1";
while (<IN1>){
	chomp ();
	@line = split/\s+/,$_; 
	$window[$line[1]][0]=$line[0];
	$window[$line[1]][1]=$line[1]; 
	$window[$line[1]][2]=$line[2];
	$type=0;
	$type_num[$line[1]][1]=0;$type_num[$line[1]][2]=0;$type_num[$line[1]][3]=0;$type_num[$line[1]][4]=0;
	foreach $i(3..$#line){
		if ($line[$i]=~/(\w+)\:(\d+)/){
			$type=$2;
			$window[$line[1]][2*$i-3]=$1;
			$window[$line[1]][2*$i-2]=$2;
			$type_num[$line[1]][$type]++;
			$sample_name[$line[1]][$type][$type_num[$line[1]][$type]]=$window[$line[1]][2*$i-3];
			$key="$line[1]	$window[$line[1]][2*$i-3]";
			$sample_type{$key}=$window[$line[1]][2*$i-2]; 
		}
	}
	if($type==3 or $type==4){
		push @more_marker_hang,$line[1];
		push @more_type,$type;
	}else{
		push @less_marker_hang,$line[1];
		push @less_type,$type;
	}
	push @all_marker_hang,$line[1];
	push @all_type,$type;
}
####选出可以聚类成3个或者4个group的位置进行分类统一###
foreach $i(1..$#more_marker_hang){
	@mapping_type=();
	foreach $j(1..$more_type[$i]){
		$missmap_type=0;@map_num=();
		foreach $m(1..$more_type[$i-1]){
			$map_num=0;
			foreach $n(1..$type_num[$more_marker_hang[$i]][$j]){
				foreach $o(1..$type_num[$more_marker_hang[$i-1]][$m]){
					if($sample_name[$more_marker_hang[$i]][$j][$n] eq $sample_name[$more_marker_hang[$i-1]][$m][$o]){
						$map_num++;
					}
				}
			}push @map_num,$map_num;
			if ($map_num/min($type_num[$more_marker_hang[$i]][$j],$type_num[$more_marker_hang[$i-1]][$m])>=0.6){
				$value_key="";
				foreach $o(1..$type_num[$more_marker_hang[$i-1]][$m]){
					$ii=$more_marker_hang[$i-1];
					$key="$ii	$sample_name[$more_marker_hang[$i-1]][$m][$o]";
					if (defined($sample_type{$key})){
						$value_key=$sample_type{$key};
						last;
					}
				}
				foreach $n(1..$type_num[$more_marker_hang[$i]][$j]){
					$key="$more_marker_hang[$i]	$sample_name[$more_marker_hang[$i]][$j][$n]";
					$sample_type{$key}="$value_key";
				}
				if(grep { $_ eq $value_key } @mapping_type){
				}else{
					push @mapping_type,$value_key;
				}
			}else{
				$missmap_type++;$value_key="";$value_key1="";$value_key2="";$value_key3="";
				if ($missmap_type==3 and $more_type[$i-1]==3){
					if(max(@map_num)<=10){
						$ii=$more_marker_hang[$i]-1;
						foreach $o(1..$type_num[$more_marker_hang[$i-1]][$m]){
							$ii=$more_marker_hang[$i-1];
							$key="$ii	$sample_name[$more_marker_hang[$i-1]][1][$o]";
							if (defined($sample_type{$key})){
								$value_key1=$sample_type{$key};
								last;
							}
						}
						foreach $o(1..$type_num[$more_marker_hang[$i-1]][$m]){
							$ii=$more_marker_hang[$i-1];
							$key="$ii	$sample_name[$more_marker_hang[$i-1]][2][$o]";
							if (defined($sample_type{$key})){
								$value_key2=$sample_type{$key};
								last;
							}
						}
						foreach $o(1..$type_num[$more_marker_hang[$i-1]][$m]){
							$ii=$more_marker_hang[$i-1];
							$key="$ii	$sample_name[$more_marker_hang[$i-1]][3][$o]";
							if (defined($sample_type{$key})){
								$value_key3=$sample_type{$key};
								last;
							}
						}
						@arr=($value_key1,$value_key2,$value_key3);
						foreach $arr(1..4){
							if(grep { $_ eq $arr } @arr){
							}else{
								$value_key=$arr;
							}
						}
						foreach $n(1..$type_num[$more_marker_hang[$i]][$j]){
							$key="$more_marker_hang[$i]	$sample_name[$more_marker_hang[$i]][$j][$n]";
							$sample_type{$key}="$value_key";
						}
						if(grep { $_ eq $value_key } @mapping_type){
						}else{
							push @mapping_type,$value_key;
						}	
					}
				}
			}
		}
	}
	
	if($#mapping_type==$more_type[$i]-1){
	}else{
		$geno{$more_marker_hang[$i]}="False";
		foreach $j(1..$more_type[$i]){
			foreach $n(1..$type_num[$more_marker_hang[$i]][$j]){
				$key="$more_marker_hang[$i]	$sample_name[$more_marker_hang[$i]][$j][$n]";
				$sample_type{$key}=$j;
			}
		}
	}
}
foreach $i(0..$#more_marker_hang){
	if(defined $geno{$more_marker_hang[$i]}){
		foreach $j(1..$more_type[$i]){
			foreach $n(1..$type_num[$more_marker_hang[$i]][$j]){
				$key="$more_marker_hang[$i]	$sample_name[$more_marker_hang[$i]][$j][$n]";
				$sample_type{$key}="";
			}
		}
	}else{
		push @select_more_marker_hang,$more_marker_hang[$i];
		push @select_more_type,$more_type[$i];
		
	}
}

					
###根据之前统一聚类的结果对聚成一类或者两类的位子进行统一###
foreach $i(0..$#less_marker_hang){
	@mapping_type=();
	foreach $j(0..$#select_more_marker_hang){
		if(max($window[$less_marker_hang[$i]][2],$window[$select_more_marker_hang[$j]][2])-min($window[$less_marker_hang[$i]][2],$window[$select_more_marker_hang[$j]][2])<=3000000){
			foreach $m(1..$less_type[$i]){
				$missmap_type=0;@map_num=();
				foreach $n(1..$select_more_type[$j]){
					$map_num=0;
					foreach $o(1..$type_num[$less_marker_hang[$i]][$m]){
						foreach $p(1..$type_num[$select_more_marker_hang[$j]][$n]){
							if($sample_name[$less_marker_hang[$i]][$m][$o] eq $sample_name[$select_more_marker_hang[$j]][$n][$p]){
								$map_num++;
							}
						}
					} 
					#print "$less_marker_hang[$i]	$select_more_marker_hang[$j]	$m	$n	$map_num	$type_num[$less_marker_hang[$i]][$m]	$type_num[$select_more_marker_hang[$j]][$n]\n";
					push @map_num,$map_num;
					if ($map_num/min($type_num[$less_marker_hang[$i]][$m],$type_num[$select_more_marker_hang[$j]][$n])>=0.9){
						$value_key="";
						foreach $p(1..$type_num[$select_more_marker_hang[$j]][$n]){
							$key="$select_more_marker_hang[$j]	$sample_name[$select_more_marker_hang[$j]][$n][$p]";
							if (defined($sample_type{$key})){
								$value_key=$sample_type{$key};
								last;
							}
						}
						foreach $o(1..$type_num[$less_marker_hang[$i]][$m]){
							$key="$less_marker_hang[$i]	$sample_name[$less_marker_hang[$i]][$m][$o]";
							$sample_type{$key}="$value_key";
						}
						if(grep { $_ eq $value_key } @mapping_type){
						}else{
						push @mapping_type,$value_key;
						}
					}else{
						$missmap_type++;$value_key="";$value_key1="";$value_key2="";$value_key3="";
						if($missmap_type==3 and $select_more_type[$j]==3){
							if(max(@map_num)<=10){
								foreach $p(1..$type_num[$select_more_marker_hang[$j]][$n]){
									$key="$select_more_marker_hang[$j]	$sample_name[$select_more_marker_hang[$j]][1][$p]";
									if (defined($sample_type{$key})){
										$value_key1=$sample_type{$key};
										last;
									}
								}
								foreach $p(1..$type_num[$select_more_marker_hang[$j]][$n]){
									$key="$select_more_marker_hang[$j]	$sample_name[$select_more_marker_hang[$j]][2][$p]";
									if (defined($sample_type{$key})){
										$value_key1=$sample_type{$key};
										last;
									}
								}
								foreach $p(1..$type_num[$select_more_marker_hang[$j]][$n]){
									$key="$select_more_marker_hang[$j]	$sample_name[$select_more_marker_hang[$j]][3][$p]";
									if (defined($sample_type{$key})){
										$value_key1=$sample_type{$key};
										last;
									}
								}
								@arr=($value_key1,$value_key2,$value_key3);
								foreach $arr(1..4){
									if(grep { $_ eq $arr } @arr){
									}else{
										$value_key=$arr;
									}
								}
								foreach $o(1..$type_num[$less_marker_hang[$i]][$m]){
									$key="$less_marker_hang[$i]	$sample_name[$less_marker_hang[$i]][$m][$o]";
									$sample_type{$key}="$value_key";
								}						
								if(grep { $_ eq $value_key } @mapping_type){
								}else{
									push @mapping_type,$value_key;
								}	
							}
						}
					}
				}
			}
		}
		if($#mapping_type == $less_type[$i]-1){####当该位置所有group都已经统一跳出这个循环，进入下一行的匹配
			last;
		}elsif($window[$select_more_marker_hang[$j]][2]-$window[$less_marker_hang[$i]][2]>3000000){
			foreach $m(1..$less_type[$i]){
				foreach $o(1..$type_num[$less_marker_hang[$i]][$m]){
					$key="$less_marker_hang[$i]	$sample_name[$less_marker_hang[$i]][$m][$o]";
					$sample_type{$key}="";
				}
			}
			last;
		}elsif($window[$less_marker_hang[$i]][2]<$window[$select_more_marker_hang[0]][2]){
			foreach $m(1..$less_type[$i]){
				foreach $o(1..$type_num[$less_marker_hang[$i]][$m]){
					$key="$less_marker_hang[$i]	$sample_name[$less_marker_hang[$i]][$m][$o]";
					$sample_type{$key}="";
				}
			}
			last;
		}elsif($window[$less_marker_hang[$i]][2]>$window[$select_more_marker_hang[$#select_more_marker_hang]][2]){
			foreach $m(1..$less_type[$i]){
				foreach $o(1..$type_num[$less_marker_hang[$i]][$m]){
					$key="$less_marker_hang[$i]	$sample_name[$less_marker_hang[$i]][$m][$o]";
					$sample_type{$key}="";
				}
			}
			last;
		}
		
	}
}
open OUT,">unification-Chr$chr.genotype";
foreach $i(0..$#all_marker_hang){
	print OUT "$all_marker_hang[$i]\t$window[$all_marker_hang[$i]][2]";
	foreach $j(@array){
		$key="$all_marker_hang[$i]	$j";
		if(defined($sample_type{$key})){
			print OUT "\t$j:$sample_type{$key}";
		}else{
			print OUT "\t$j:";
		}
	}
	print OUT "\n";
}
close OUT;
close IN;
close IN1;
open IN2,"unification-Chr$chr.genotype" or die "can't open the input 2 file";
while (<IN2>){
	chomp();
	@marker = split/\s+/,$_; 
	$window[$marker[0]][0]=$marker[0];
	$window[$marker[0]][1]=$marker[1];
	foreach $i(2..$#marker){
		@sample=();
		@sample=split/\:/,$marker[$i];
		$sample_name[$marker[0]][$i]=$sample[0];
		if ($sample[1]=~/[1 2 3 4]/){
			$window[$marker[0]][$i]=$sample[1];
		}else{
			$window[$marker[0]][$i]="0";
		}
	}
}
foreach $i(2..$#marker){ 
	@marker_num=();@marker_type=();
	foreach $j(1..$marker[0]){
		if ($window[$j][$i]!=0){
			push @marker_num,$window[$j][0];
			push @marker_type,$window[$j][$i];
		}
	}
	foreach $m(0..$#marker_num){
		if ($m==0 and $marker_num[0]>1){
			foreach $n(1..$marker_num[0]-1){
				$window[$n][$i]="$marker_type[0]";
			}
		}else{
			if($marker_num[$m]-$marker_num[$m-1]>1){
				if ($marker_type[$m] == $marker_type[$m-1] ){
					$star1=$marker_num[$m-1]+1;$end1=$marker_num[$m]-1;
					foreach $n($star1..$end1){
						$window[$n][$i]=$marker_type[$m-1];
					}
				}else{
					$star1=$marker_num[$m-1]+1;
					$end1=int(($marker_num[$m]-$marker_num[$m-1])/2)+$marker_num[$m-1];
					$star2=int(($marker_num[$m]-$marker_num[$m-1])/2)+$marker_num[$m-1]+1;
					$end2=$marker_num[$m]-1;
					foreach $n($star1..$end1){
						$window[$n][$i]=$marker_type[$m-1];
					}
					foreach $n($star2..$end2){
						$window[$n][$i]=$marker_type[$m];
					}
				}
			}
		}
	}
	if($marker_num[$#marker_num] != $marker[0]){
		foreach $n($marker_num[$#marker_num]+1..$marker[0]){
			$window[$n][$i]=$marker_type[$#marker_num];
		} 
	}
	
}
open OUT2,">Imputation-Chr$chr.genotype";
foreach $j(1..$marker[0]){
	print OUT2 "$window[$j][0]\t$window[$j][1]";
	foreach $i(2..$#marker){
		if($window[$j][$i]==1){
			print OUT2 "\tAC";
		}elsif($window[$j][$i]==2){
			print OUT2 "\tAD";
		}elsif($window[$j][$i]==3){
			print OUT2 "\tBC";
		}elsif($window[$j][$i]==4){
			print OUT2 "\tBD";
		}else{
			print OUT2 "\t$window[$j][$i]";
		}
		
	}
	print OUT2 "\n";
}
close IN2;
close OUT2;
