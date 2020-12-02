$chr=$ARGV[1];$cuttree= $ARGV[3];
$marker_num=$ARGV[2];$min_window_marker_num=$ARGV[4];$max_window_marker_num=$ARGV[5];
use List::Util qw/max min/;
use Statistics::Basic qw(:all);
$hang=0;$window_num=1;$window_hang=1;
open IN1,"$ARGV[0].geno" or die "can't open the input file";
while (<IN1>){
	chomp();
	$hang++;
	$line = $_;
	if ($hang>1){
		if ($window_hang<=$marker_num){
			$key="$window_num-$window_hang";
			$marker{$key}=$line;
			$window_hang++;
		}else{
			$window_hang=1;
			$window_num++;
			$key="$window_num-$window_hang";
			$marker{$key}=$line;
			$window_hang++;
		}
	}
}
$last_window_hang=$window_hang-1;
foreach $i(1..$window_num){
	if ($i<$window_num){
		foreach $j(1..$marker_num){
			$key="$i-$j";
		}
	}elsif ($i==$window_num){
		if ($last_window_hang<$marker_num*0.6){
			foreach $j($marker_num+1..$last_window_hang+$marker_num){
				$ii=$i-1;
				$key1="$ii".'-'."$j";
				$jj=$j-$marker_num;
				$key2="$i-"."$jj";
				$marker{$key1}=$marker{$key2};
			}
			$window_num=$window_num-1;
			$last_window_hang=$last_window_hang+$marker_num;
		}else{
			$window_num=$window_num;
			$last_window_hang=$last_window_hang;
		}
	}
}

open OUTFILE,">R2-$ARGV[0].txt";
foreach $i(1..$window_num){
	@snp = ();
	if ($i<$window_num){
		$window_marker_num=$marker_num;
	}elsif($i==$window_num){
		$window_marker_num=$last_window_hang;
	}
	foreach $j(1..$window_marker_num){
		$key="$i-$j";
		$hang_marker = $marker{$key};
		if($hang_marker=~/^Chr(\d+)\s+(\d+)/){
			$chr = $1;$pos = $2;$key1="Chr$chr".'_'."$pos";
			$hang_snp{$key1}=$hang_marker;
			@array=split/\s+/,$hang_marker;
			foreach $ii(0..5){
				$snp[$j][$ii] = $array[$ii];
			}
			foreach $ii(8..$#array){
				if($array[$ii] eq "H"){
					$snp[$j][$ii] = 1; 
				}elsif($array[$ii] eq "-"){
					$snp[$j][$ii] = undef; 
				}else{
					$snp[$j][$ii] = 0;
				}
			}
		}
		
	}
	$out="$ARGV[0]-"."$i".'-r2';
	open OUT,">$out.txt" or "die can't open the out file";
	foreach (1..$window_marker_num){
		if ($_==1){
			print OUT "$snp[$_][0]_$snp[$_][1]";
		}else{
			print OUT "\t$snp[$_][0]_$snp[$_][1]";
		}
	}
	print OUT "\n";
	foreach $ii(1..$window_marker_num){
		$num_correlation =0;@id_correlation=();
		foreach $jj(1..$window_marker_num){
			if ($jj==1){
				if ($ii==$jj){
					$correlation[$ii][$jj] =sprintf "%0.4f",0;
					print OUT "$correlation[$ii][$jj]";
				}elsif($ii<$jj){
					@v1 =();@v2 = ();
					foreach (8..$#array){
						push @v1,$snp[$ii][$_];
						push @v2,$snp[$jj][$_];
					}
					$v1 = vector @v1;
					$v2 = vector @v2;
					($f1, $f2) = handle_missing_values($v1, $v2);
					$correlation=correlation ($f1, $f2);
					$correlation[$ii][$jj]=sprintf "%0.4f",(1-$correlation*$correlation) ;
					print OUT "$correlation[$ii][$jj]";
				}elsif($ii>$jj){
					$correlation[$ii][$jj]=$correlation[$jj][$ii];
					print OUT "$correlation[$ii][$jj]";
				}
			}else{
				if ($ii==$jj){
					$correlation[$ii][$jj] =sprintf "%0.4f",0;
					print OUT "\t$correlation[$ii][$jj]";
				}elsif($ii<$jj){
					@v1 =();@v2 = ();
					foreach (8..$#array){
						push @v1,$snp[$ii][$_];
						push @v2,$snp[$jj][$_];
					}
					$v1 = vector @v1;
					$v2 = vector @v2;
					($f1, $f2) = handle_missing_values($v1, $v2);
					$correlation=correlation ($f1, $f2);
					$correlation[$ii][$jj]=sprintf "%0.4f",(1-$correlation*$correlation) ;
					print OUT "\t$correlation[$ii][$jj]";
				}elsif($ii>$jj){
					$correlation[$ii][$jj]=$correlation[$jj][$ii];
					print OUT "\t$correlation[$ii][$jj]";
				}
			}
		}
			print OUT "\n";
	}
	$temp='Rscript cluster.R '."$out"." $cuttree";
	print "$temp\n";
	system ("$temp");
	$cluster_file = "$out"."-$cuttree";
	@cluster_type=();%cluster=();@cluster_name=();@select_type=();
	open IN3,"$cluster_file.txt";
	while(<IN3>){
		chomp();	
		@cluster_snp=split/\s+/,$_;
		$cluster{$cluster_snp[0]} = $cluster_snp[1];
		push @cluster_type,$cluster_snp[1];
		push @cluster_name,$cluster_snp[0];
	}
	@cluster_type = sort {$a <=> $b} @cluster_type;
	$star =0;
	foreach (@cluster_type){
		if ($_ == $star){
			$star_num++;
		}else{
			if ($star_num >= $min_window_marker_num and $star_num <= $max_window_marker_num ){
				push @select_type,$star;
			}
			$star++;
			$star_num=1;
		}
	}
	foreach $i(@select_type){
		foreach $j(@cluster_name){
			if ($i==$cluster{$j}){
				print OUTFILE "$hang_snp{$j}\n";
			}
		}
		print OUTFILE "\n";
	}
}
close IN1;
close IN3;
close OUT;
close OUTFILE;	
		
		
		
		
		
		
		
		
		














