my (@line,$chr,$win_length);
$win_length=$ARGV[0];
######################
$min_snp_num = $ARGV[1];
$chromosome = $ARGV[3];$window = 1;$all_hang_num= 0;$star_hang_num{1}= 2;
open IN2,"$ARGV[2]";
while(<IN2>){
	chomp();
	$hang=$_;
	if($hang=~/^Chr(\d+)\s+(\d+)/){
		$chr =$1;$pos=$2;
		if ($chr == $chromosome){
			$window_num = int($pos/$win_length)+1;
			if ($window_num == $window ){
				$window_snp++;$all_hang_num++;
			}else{
				$all_hang_num++;
				$window_snp_num{$window} = $window_snp;  
				$window_snp = 1;$window++;
			}
		}
	}elsif($_=~/^#CHROM/){
		@id = split/\s+/,$_;
		shift (@id);shift (@id);shift (@id);shift (@id);
	}
}
$window_snp_num{$window} = $window_snp; 
#######################
$c = 0; $combin_num=0;
foreach $i(1..$window){
	$c = $window_snp_num{$i}+$c;
	if ($c < $min_snp_num ){
		$combin_num++;
	}else{
		$chr_window++;
		$star{$chr_window}=($i-1-$combin_num)*$win_length;
		$end{$chr_window}=$i*$win_length-1;
		$chr_window_snp_num{$chr_window}=$c;
		$combin_num=0;$c=0;
	}
}
if ($c<$min_snp_num){
	$chr_window_snp_num{$chr_window}=$c+$chr_window_snp_num{$chr_window};
	$end{$chr_window}=$window*$win_length-1;
}else{
	$chr_window++;
	$chr_window_snp_num{$chr_window}=$c;
	$star{$chr_window}=($window-1-$combin_num)*$win_length;
	$end{$chr_window}=$window*$win_length-1;
}
###########################	
$chr_window_num=1;
open IN2,"$ARGV[2]";
while(<IN2>){
	chomp();
	$array=$_;
	if($array=~/^Chr(\d+)\s+(\d+)/){
		$chr = $1;$pos = $2;
		@array=split/\s+/,$array;
		if ($pos>=$star{$chr_window_num} and $pos<=$end{$chr_window_num} ){ 
			$a++;
			if($a<$chr_window_snp_num{$chr_window_num}){
				foreach $i(4..$#array){ 
					if ($array[$i] eq $array[2]){
						$snp[$a][$i-4] = 0;  
					}elsif($array[$i] eq $array[3]){
						$snp[$a][$i-4] = 2;   
					}elsif($array[$i] eq "H"){
						$snp[$a][$i-4] = 1; 
					}elsif($array[$i] eq "-"){
						$snp[$a][$i-4] = 9; 
					}
				}
			}elsif($a==$chr_window_snp_num{$chr_window_num}){ 
				foreach $i(4..$#array){ 
					if ($array[$i] eq $array[2]){
						$snp[$a][$i-4] = 0;  
					}elsif($array[$i] eq $array[3]){
						$snp[$a][$i-4] = 2;   
					}elsif($array[$i] eq "H"){
						$snp[$a][$i-4] = 1; 
					}elsif($array[$i] eq "-"){
						$snp[$a][$i-4] = 9; 
					}
				}
				print "Calculating Chr$chromosome $chr_window_num kinship\n";
				$filename = 'Chr'."$chromosome".'_'."$chr_window_num".'_'."$end{$chr_window_num}"."-$ARGV[4].kinf";
				open OUT,">$filename";
				foreach $i(0..$#id){
					if ($i==0){
						print OUT "$id[$i]";
					}else{
						print OUT "\t$id[$i]";
					}
				}
				print OUT "\n";
				foreach $i(0..$#array-4){ 
					foreach $j(0..$#array-4){
						$Same=0;$Total=0;
						if ($j == 0 ){
							if($i == $j){
								$Kin[$i][$j]=sprintf "%.4f",0;print OUT "$Kin[$i][$j]";
							}elsif($i > $j){
								$Kin[$i][$j]=$Kin[$j][$i];
								print OUT "$Kin[$i][$j]";
							}elsif($i < $j){
								foreach $n(1..$a){
									$Geno_i=$snp[$n][$i];
									$Geno_j=$snp[$n][$j];
									$Total=$Total+1;
									if ($Geno_i == 9 or $Geno_j == 9){	
										$Total=$Total-1;
									}elsif ($Geno_i == $Geno_j){
										$Same=$Same+1;
									}elsif($Geno_i == 1 or  $Geno_j == 1){
										$Same=$Same+0.5;
									}
								}
								$Kin[$i][$j]=sprintf "%.4f",1-$Same/$Total; 
								print OUT "$Kin[$i][$j]";
							}
						}else{
							if($i == $j){
								$Kin[$i][$j]=sprintf "%.4f",0;print OUT "\t$Kin[$i][$j]";
							}elsif($i > $j){
								$Kin[$i][$j]=$Kin[$j][$i];
								print OUT "\t$Kin[$i][$j]";
							}elsif($i < $j){
								foreach $n(1..$a){
									$Geno_i=$snp[$n][$i];
									$Geno_j=$snp[$n][$j];
									$Total=$Total+1;
									if ($Geno_i == 9 or $Geno_j == 9){	
										$Total=$Total-1;
									}elsif ($Geno_i == $Geno_j){
										$Same=$Same+1;
									}elsif($Geno_i == 1 or  $Geno_j == 1){
										$Same=$Same+0.5;
									}
								}
								if ($Total >= 1){
									$Kin[$i][$j]=sprintf "%.4f",1-$Same/$Total; 
								 }elsif($Total == 0){
									$Kin[$i][$j]=1.0000;
								}
								print OUT "\t$Kin[$i][$j]";
							}
						}
					}
					print OUT "\n";
				}
				close OUT;
				$a=0;@snp = ();@Kin=();$chr_window_num++;
			}
		}
	}
}
close IN2;
