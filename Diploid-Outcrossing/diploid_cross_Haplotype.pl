open OUT,">$ARGV[1]_Haplotype.geno";
open IN1,"$ARGV[0]" or die "can't open the input file 1";
$chromosome=1;
while (<IN1>){
	chomp();
	if ($.==1){
		@sample=split/\t/,$_;
		 shift(@sample);shift(@sample);
	}else{
		$hang++;
		$window_num++;
		@window=split/\t/,$_;
		$geno{$window_num}=$_;
		$chr{$window_num}=$window[0];
		$pos{$window_num}=$window[1];
		if ($window[0]==$chromosome){
		}else{
			$hang_num{$chromosome}=$hang;
			$chromosome++;
		}
			
	}
}
$hang_num{$chromosome}=$hang;
$chromosome=1;$pos{0}=0;
foreach $i(1..$window_num){
	if ($i<$hang_num{$chromosome}){
		if($i==1){
			$pos_star{$i}=0;
			$pos_end{$i}=$pos{$i}*2000000-1;
		}else{
			$pos_star{$i}=$pos_end{$i-1}+1;
			$pos_end{$i}=($pos{$i}*2-($pos_star{$i}+1)/1000000)*1000000;
		}
	}else{
		$chromosome++;
		$pos_star{$i}=0;
		$pos_end{$i}=$pos{$i}*2000000-1;
	}
}


$i=1;
$m=0;

open IN2,"$ARGV[1].vcf" or die "can't open the input file 2";
while (<IN2>){
	chomp();
	@snp=split/\s+/,$_;
	if ($snp[0] eq "#CHROM"){
		foreach $j(9..$#snp){   
			push @vcf_sample,$snp[$j];
			if ($snp[$j] eq $ARGV[2]){
				$parent1=$j;
			}elsif($snp[$j] eq $ARGV[3]){
				$parent2=$j;
			}
		}
	}else{
		if ($snp[0]=~/^Chr(\d+)/){
			$chr =$1;
			if ($chr{$i} == $chr and $snp[1] >= $pos_star{$i} and $snp[1] <= $pos_end{$i}){$marker_ture++;
				if($marker_ture==1){
					$star=$i;
				}
				if (length($snp[3])==1 and length ($snp[4])==1 and $snp[5]>=300 and $snp[6] ne "LowQual" and $snp[$parent1] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/ and $snp[$parent2] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
					$m++;
					if ($m==1){
						@geno="";$geno1="";$geno2="";$geno3="";$geno4="";
						@geno1="";@geno2="";@geno3="";@geno4="";
						@geno1_num="";@geno2_num="";@geno3_num="";@geno4_num="";
						@geno=split/\s+/,$geno{$i};
						shift(@geno);shift(@geno);
						foreach $j(0..$#geno){
							if($geno[$j] eq "AC"){
								$geno1 .= "$sample[$j] ";
							}elsif($geno[$j] eq "AD"){
								$geno2 .= "$sample[$j] ";
							}elsif($geno[$j] eq "BC"){
								$geno3 .= "$sample[$j] ";
							}elsif($geno[$j] eq "BD"){
								$geno4 .= "$sample[$j] ";
							}
						}
						@geno1=split/\s+/,$geno1;@geno2=split/\s+/,$geno2;
						@geno3=split/\s+/,$geno3;@geno4=split/\s+/,$geno4;
						foreach $o(@geno1){
							foreach $p(0..$#vcf_sample){
								if ($o eq $vcf_sample[$p]){
									push @geno1_num,$p+9;
								}
							}
						}
						foreach $o(@geno2){
							foreach $p(0..$#vcf_sample){
								if ($o eq $vcf_sample[$p]){
									push @geno2_num,$p+9;
								}
							}
						}
						foreach $o(@geno3){
							foreach $p(0..$#vcf_sample){
								if ($o eq $vcf_sample[$p]){
									push @geno3_num,$p+9;
								}
							}
						}
						foreach $o(@geno4){
							foreach $p(0..$#vcf_sample){
								if ($o eq $vcf_sample[$p]){
									push @geno4_num,$p+9;
								}
							}
						}
						$ref_1=0;$alt_1=0;$ref_2=0;$alt_2=0;$ref_3=0;$alt_3=0;$ref_4=0;$alt_4=0;
						foreach (@geno1_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_1=$ref_1+$3;
								$alt_1=$alt_1+$4;
							}
						}
						foreach (@geno2_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_2=$ref_2+$3;
								$alt_2=$alt_2+$4;
							}
						}foreach (@geno3_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_3=$ref_3+$3;
								$alt_3=$alt_3+$4;
							}
						}foreach (@geno4_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_4=$ref_4+$3;
								$alt_4=$alt_4+$4;
							}
						}
						if ($snp[$parent1] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/) {
							$parent1_num1=$3;$parent1_num2=$4;
							if ($snp[$parent2] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/ ){
								$marker[$i][$m][0]="$snp[0]";$marker[$i][$m][1]="$snp[1]";$marker[$i][$m][2]="$snp[3]";$marker[$i][$m][3]="$snp[4]";
								$marker[$i][$m][4]=$parent1_num1;$marker[$i][$m][5]=$parent1_num2;$marker[$i][$m][6]=$3;$marker[$i][$m][7]=$4;
								$marker[$i][$m][8]=$ref_1;$marker[$i][$m][9]=$alt_1;$marker[$i][$m][10]=$ref_2;$marker[$i][$m][11]=$alt_2;
								$marker[$i][$m][12]=$ref_3;$marker[$i][$m][13]=$alt_3;$marker[$i][$m][14]=$ref_4;$marker[$i][$m][15]=$alt_4;
							}
						}
					}else{
						$ref_1=0;$alt_1=0;$ref_2=0;$alt_2=0;$ref_3=0;$alt_3=0;$ref_4=0;$alt_4=0;
						foreach (@geno1_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_1=$ref_1+$3;
								$alt_1=$alt_1+$4;
							}
						}
						foreach (@geno2_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_2=$ref_2+$3;
								$alt_2=$alt_2+$4;
							}
						}foreach (@geno3_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_3=$ref_3+$3;
								$alt_3=$alt_3+$4;
							}
						}foreach (@geno4_num){
							if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
								$ref_4=$ref_4+$3;
								$alt_4=$alt_4+$4;
							}
						}
						if ($snp[$parent1] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/) {
							$parent1_num1=$3;$parent1_num2=$4;
							if ($snp[$parent2] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/ ){
								$marker[$i][$m][0]="$snp[0]";$marker[$i][$m][1]="$snp[1]";$marker[$i][$m][2]="$snp[3]";$marker[$i][$m][3]="$snp[4]";
								$marker[$i][$m][4]=$parent1_num1;$marker[$i][$m][5]=$parent1_num2;$marker[$i][$m][6]=$3;$marker[$i][$m][7]=$4;
								$marker[$i][$m][8]=$ref_1;$marker[$i][$m][9]=$alt_1;$marker[$i][$m][10]=$ref_2;$marker[$i][$m][11]=$alt_2;
								$marker[$i][$m][12]=$ref_3;$marker[$i][$m][13]=$alt_3;$marker[$i][$m][14]=$ref_4;$marker[$i][$m][15]=$alt_4;
							}
						}
					}
				}
			}elsif ($chr{$i} == $chr and $snp[1] >= $pos_end{$i} and length($snp[3])==1 and length ($snp[4])==1 and $snp[5]>=300 and $snp[6] ne "LowQual" and $snp[$parent1] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/ and $snp[$parent2] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
				push @m,$m;
				$m=1;
				$i++;
				@geno="";$geno1="";$geno2="";$geno3="";$geno4="";
				@geno1="";@geno2="";@geno3="";@geno4="";
				@geno1_num="";@geno2_num="";@geno3_num="";@geno4_num="";
				@geno=split/\s+/,$geno{$i};
				shift(@geno);shift(@geno);
				foreach $j(0..$#geno){
					if($geno[$j] eq "AC"){
						$geno1 .= "$sample[$j] ";
					}elsif($geno[$j] eq "AD"){
						$geno2 .= "$sample[$j] ";
					}elsif($geno[$j] eq "BC"){
						$geno3 .= "$sample[$j] ";
					}elsif($geno[$j] eq "BD"){
						$geno4 .= "$sample[$j] ";
					}
				}
				@geno1=split/\s+/,$geno1;@geno2=split/\s+/,$geno2;
				@geno3=split/\s+/,$geno3;@geno4=split/\s+/,$geno4;
				foreach $o(@geno1){
					foreach $p(0..$#vcf_sample){
						if ($o eq $vcf_sample[$p]){
							push @geno1_num,$p+9;
						}
					}
				}
				foreach $o(@geno2){
					foreach $p(0..$#vcf_sample){
						if ($o eq $vcf_sample[$p]){
							push @geno2_num,$p+9;
						}
					}
				}
				foreach $o(@geno3){
					foreach $p(0..$#vcf_sample){
						if ($o eq $vcf_sample[$p]){
							push @geno3_num,$p+9;
						}
					}
				}
				foreach $o(@geno4){
					foreach $p(0..$#vcf_sample){
						if ($o eq $vcf_sample[$p]){
							push @geno4_num,$p+9;
						}
					}
				}
				$ref_1=0;$alt_1=0;$ref_2=0;$alt_2=0;$ref_3=0;$alt_3=0;$ref_4=0;$alt_4=0;
				foreach (@geno1_num){
					if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
						$ref_1=$ref_1+$3;
						$alt_1=$alt_1+$4;
					}
				}
				foreach (@geno2_num){
					if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
						$ref_2=$ref_2+$3;
						$alt_2=$alt_2+$4;
					}
				}foreach (@geno3_num){
					if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
						$ref_3=$ref_3+$3;
						$alt_3=$alt_3+$4;
					}
				}foreach (@geno4_num){
					if ($snp[$_]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
						$ref_4=$ref_4+$3;
						$alt_4=$alt_4+$4;
					}
				}
				if ($snp[$parent1] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/) {
					$parent1_num1=$3;$parent1_num2=$4;
					if ($snp[$parent2] =~/(\d)\/(\d)\:(\d+)\,(\d+)\:/ ){
						$marker[$i][$m][0]="$snp[0]";$marker[$i][$m][1]="$snp[1]";$marker[$i][$m][2]="$snp[3]";$marker[$i][$m][3]="$snp[4]";
						$marker[$i][$m][4]=$parent1_num1;$marker[$i][$m][5]=$parent1_num2;$marker[$i][$m][6]=$3;$marker[$i][$m][7]=$4;
						$marker[$i][$m][8]=$ref_1;$marker[$i][$m][9]=$alt_1;$marker[$i][$m][10]=$ref_2;$marker[$i][$m][11]=$alt_2;
						$marker[$i][$m][12]=$ref_3;$marker[$i][$m][13]=$alt_3;$marker[$i][$m][14]=$ref_4;$marker[$i][$m][15]=$alt_4;
					}
				}
			}elsif($chr{$i} != $chr){
				$i++;
			}
		}
	}
}		
push @m,$m;
foreach $j($star..$i){
	if($j>1){push @jj,$jj;}
	$jj=0;
	foreach $m(1..$m[$j-$star]){
		if (($marker[$j][$m][4]+$marker[$j][$m][5])>= 30 and ($marker[$j][$m][6]+$marker[$j][$m][7])>=30 and ($marker[$j][$m][8]+$marker[$j][$m][9])>=30  and ($marker[$j][$m][11]+$marker[$j][$m][10])>=30 and ($marker[$j][$m][12]+$marker[$j][$m][13])>=30 and ($marker[$j][$m][14]+$marker[$j][$m][15])>=30){
			if ($marker[$j][$m][4]/($marker[$j][$m][4]+$marker[$j][$m][5])>=0.9){
				$parent1_geno=$marker[$j][$m][2];
			}elsif($marker[$j][$m][5]/($marker[$j][$m][4]+$marker[$j][$m][5])>=0.9){
				$parent1_geno=$marker[$j][$m][3];
			}elsif($marker[$j][$m][5]/($marker[$j][$m][4]+$marker[$j][$m][5])>=0.3 and $marker[$j][$m][5]/($marker[$j][$m][4]+$marker[$j][$m][5])<=0.7){
				$parent1_geno="H";
			}else{
				$parent1_geno="none";
			}
			if ($marker[$j][$m][6]/($marker[$j][$m][6]+$marker[$j][$m][7])>=0.9){
				$parent2_geno=$marker[$j][$m][2];
			}elsif($marker[$j][$m][7]/($marker[$j][$m][6]+$marker[$j][$m][7])>=0.9){
				$parent2_geno=$marker[$j][$m][3];
			}elsif($marker[$j][$m][6]/($marker[$j][$m][6]+$marker[$j][$m][7])>=0.3 and $marker[$j][$m][6]/($marker[$j][$m][6]+$marker[$j][$m][7])<=0.7){
				$parent2_geno="H";
			}else{
				$parent2_geno="none";
			}
			if ($marker[$j][$m][8]/($marker[$j][$m][8]+$marker[$j][$m][9])>=0.9){
				$genotype1=$marker[$j][$m][2];
			}elsif($marker[$j][$m][9]/($marker[$j][$m][8]+$marker[$j][$m][9])>=0.9){
				$genotype1=$marker[$j][$m][3];
			}elsif($marker[$j][$m][8]/($marker[$j][$m][8]+$marker[$j][$m][9])>=0.3 and $marker[$j][$m][8]/($marker[$j][$m][8]+$marker[$j][$m][9])<=0.7){
				$genotype1="H";
			}else{
				$genotype1="none";
			}
			if ($marker[$j][$m][10]/($marker[$j][$m][10]+$marker[$j][$m][11])>=0.9){
				$genotype2=$marker[$j][$m][2];
			}elsif($marker[$j][$m][11]/($marker[$j][$m][10]+$marker[$j][$m][11])>=0.9){
				$genotype2=$marker[$j][$m][3];
			}elsif($marker[$j][$m][10]/($marker[$j][$m][10]+$marker[$j][$m][11])>=0.3 and $marker[$j][$m][10]/($marker[$j][$m][10]+$marker[$j][$m][11])<=0.7){
				$genotype2="H";
			}else{
				$genotype2="none";
			}
			if ($marker[$j][$m][12]/($marker[$j][$m][12]+$marker[$j][$m][13])>=0.9){
				$genotype3=$marker[$j][$m][2];
			}elsif($marker[$j][$m][13]/($marker[$j][$m][12]+$marker[$j][$m][13])>=0.9){
				$genotype3=$marker[$j][$m][3];
			}elsif($marker[$j][$m][12]/($marker[$j][$m][12]+$marker[$j][$m][13])>=0.3 and $marker[$j][$m][12]/($marker[$j][$m][12]+$marker[$j][$m][13])<=0.7){
				$genotype3="H";
			}else{
				$genotype3="none";
			}
			if ($marker[$j][$m][14]/($marker[$j][$m][14]+$marker[$j][$m][15])>=0.9){
				$genotype4=$marker[$j][$m][2];
			}elsif($marker[$j][$m][15]/($marker[$j][$m][14]+$marker[$j][$m][15])>=0.9){
				$genotype4=$marker[$j][$m][3];
			}elsif($marker[$j][$m][14]/($marker[$j][$m][14]+$marker[$j][$m][15])>=0.3 and $marker[$j][$m][14]/($marker[$j][$m][14]+$marker[$j][$m][15])<=0.7){
				$genotype4="H";
			}else{
				$genotype4="none";
			}
			@F_geno= "";
			@F_geno=($genotype1,$genotype2,$genotype3,$genotype4);
			$heter=0;$ref=0;$alt=0;
			foreach (@F_geno){
				if ($_ eq "H"){
					$heter++;
				}elsif($_ eq $marker[$j][$m][2]){
					$ref++;
				}elsif($_ eq $marker[$j][$m][3]){
					$alt++;
				}
			}
			
			
			if($parent1_geno ne "none" and $parent2_geno ne "none" and $genotype1 ne "none"  and $genotype2 ne "none" and $genotype3 ne "none" and $genotype4 ne "none" ){
				if($parent1_geno eq $parent2_geno and $parent1_geno eq  "H")  {
					if ($ref ==1 and $alt ==1 and $heter ==2){ 
						$jj++;
						$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
						$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
						$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
					}
				}elsif($parent1_geno eq "H" and $parent2_geno eq "$marker[$j][$m][2]"){
					if($ref ==2 and $heter == 2){
						$jj++;
						$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
						$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
						$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
					}
				}elsif($parent1_geno eq "H" and $parent2_geno eq "$marker[$j][$m][3]"){
					if($ref ==2 and $heter == 2){
						$jj++;
						$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
						$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
						$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
					}
				}elsif($parent2_geno eq "H" and $parent1_geno eq "$marker[$j][$m][3]"){
					if($ref ==2 and $heter == 2){
						$jj++;
						$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
						$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
						$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
					}
				}elsif($parent2_geno eq "H" and 	$parent1_geno eq "$marker[$j][$m][2]"){
					if($ref ==2 and $heter == 2){
						$jj++;
						$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
						$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
						$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
					}
				}elsif($parent2_geno eq "$marker[$j][$m][3]" and 	$parent1_geno eq "$marker[$j][$m][2]" and $heter == 4){
					$jj++;
					$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
					$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
					$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
				}elsif($parent2_geno eq "$marker[$j][$m][2]" and 	$parent1_geno eq "$marker[$j][$m][3]" and $heter == 4){
					$jj++;
					$pos[$j][$jj][0]=$marker[$j][$m][0];$pos[$j][$jj][1]=$marker[$j][$m][1];$pos[$j][$jj][2]=$marker[$j][$m][2];
					$pos[$j][$jj][3]=$marker[$j][$m][3];$pos[$j][$jj][4]=$parent1_geno;$pos[$j][$jj][5]=$parent2_geno;
					$pos[$j][$jj][6]=$genotype1;$pos[$j][$jj][7]=$genotype2;$pos[$j][$jj][8]=$genotype3;$pos[$j][$jj][9]=$genotype4;
				}
			}
		}
	}
}
push @jj,$jj;
foreach $j($star..$i){	
	$num2=0;$num1=0;
	$parent1_1="";$parent2_1="";$geno1_1="";$geno2_1="";$geno3_1="";$geno4_1="";$geno_1="";
	$parent1_2="";$parent2_2="";$geno1_2="";$geno2_2="";$geno3_2="";$geno4_2="";$geno_2="";
	foreach $jj(1..$jj[$j-$star]){
		if ($pos[$j][$jj][4] eq "H" and $pos[$j][$jj][5] ne "H" and $num1==0){
			$num1++;$geno_1=$pos[$j][$jj][5];
			$parent1_1=$pos[$j][$jj][4];$parent2_1=$pos[$j][$jj][5];$geno1_1=$pos[$j][$jj][6];$geno2_1=$pos[$j][$jj][7];$geno3_1=$pos[$j][$jj][8];$geno4_1=$pos[$j][$jj][9];
		}elsif($pos[$j][$jj][4] ne "H" and $pos[$j][$jj][5] eq "H" and $num2==0){
			$num2++;$geno_2=$pos[$j][$jj][4];
			$parent1_2=$pos[$j][$jj][4];$parent2_2=$pos[$j][$jj][5];$geno1_2=$pos[$j][$jj][6];$geno2_2=$pos[$j][$jj][7];$geno3_2=$pos[$j][$jj][8];$geno4_2=$pos[$j][$jj][9];
		}
		if (($num1+$num2)==2){
			@site1="";@site2="";
			@site1=($geno1_1,$geno2_1,$geno3_1,$geno4_1);
			@site2=($geno1_2,$geno2_2,$geno3_2,$geno4_2);
			foreach (0..$#site1){
				if ($site1[$_] eq $geno_1){
					$site_type[$_ ]=1;
				}else{
					$site_type[$_ ]=2;
				}
			}
			foreach (0..$#site2){
				if ($site2[$_] eq $geno_2){
					$site_type[$_ ].=3;
				}else{
					$site_type[$_ ].=4;
				}
			}
			foreach (0..$#site2){
				$geno_type[$j][$_+1]=$site_type[$_ ];
			}
			last;
		}
	}
}
print OUT "Chr	Pos	Ref	Alt	Parent1	Parent2	Haplotype1	Haplotype2	Haplotype3	Haplotype4\n";
foreach $j($star..$i){
	foreach $jj(1..$jj[$j-$star]){
		@F_geno= "";
		@F_geno=($pos[$j][$jj][6],$pos[$j][$jj][7],$pos[$j][$jj][8],$pos[$j][$jj][9]);
		$heter=0;$ref=0;$alt=0;$haplotype1="";$haplotype2="";	$haplotype3="";$haplotype4="";
		foreach (@F_geno){
			if ($_ eq "H"){
				$heter++;
			}elsif($_ eq $pos[$j][$jj][2]){
				$ref++;
			}elsif($_ eq $pos[$j][$jj][3]){
				$alt++;
			}
		}
		foreach (1..4){
			if ($geno_type[$j][$_]==13){
				$AC=$pos[$j][$jj][$_+5];
			}elsif($geno_type[$j][$_]==14){
				$AD=$pos[$j][$jj][$_+5];
			}elsif($geno_type[$j][$_]==23){
				$BC=$pos[$j][$jj][$_+5];
			}elsif($geno_type[$j][$_]==24){
				$BD=$pos[$j][$jj][$_+5];
			}
		}
		if ($pos[$j][$jj][4] eq $pos[$j][$jj][5] and$pos[$j][$jj][4] eq "H"){
			if ($alt == 1 and $ref ==1 and $heter == 2){	
				if ($AC	eq $pos[$j][$jj][2]){
					$haplotype1=$pos[$j][$jj][2];$haplotype3=$pos[$j][$jj][2];
					if($BD	eq $pos[$j][$jj][3]){
						$haplotype2=$pos[$j][$jj][3];$haplotype4=$pos[$j][$jj][3];
						print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
					}
				}elsif($AD	eq $pos[$j][$jj][2]){
					$haplotype1=$pos[$j][$jj][2];$haplotype4=$pos[$j][$jj][2];
					if($BC	eq $pos[$j][$jj][3]){
						$haplotype2=$pos[$j][$jj][3];$haplotype3=$pos[$j][$jj][3];
						print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
					}
				}elsif($BC	eq $pos[$j][$jj][2]){
					$haplotype2=$pos[$j][$jj][2];$haplotype3=$pos[$j][$jj][2];
					if($AD	eq $pos[$j][$jj][3]){
						$haplotype1=$pos[$j][$jj][3];$haplotype4=$pos[$j][$jj][3];
						print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
					}
				}elsif($BD	eq $pos[$j][$jj][2]){
					$haplotype2=$pos[$j][$jj][2];$haplotype4=$pos[$j][$jj][2];
					if($AC	eq $pos[$j][$jj][3]){
						$haplotype1=$pos[$j][$jj][3];$haplotype3=$pos[$j][$jj][3];
						print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
					}
				}
			}
		}elsif($pos[$j][$jj][4] eq  $pos[$j][$jj][2] and $pos[$j][$jj][5] eq  $pos[$j][$jj][3] and $heter == 4){
			print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$pos[$j][$jj][2]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][3]\n";	
		}elsif($pos[$j][$jj][4] eq  $pos[$j][$jj][3] and $pos[$j][$jj][5] eq  $pos[$j][$jj][2] and $heter == 4){
			print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$pos[$j][$jj][3]	$pos[$j][$jj][3]	$pos[$j][$jj][2]	$pos[$j][$jj][2]\n";
		}elsif($pos[$j][$jj][4] eq   $pos[$j][$jj][2] and $pos[$j][$jj][5] eq    "H"){
			if ($ref ==2 and $heter == 2){
				$haplotype1=$pos[$j][$jj][2];$haplotype2=$pos[$j][$jj][2];
				if($AC	eq $pos[$j][$jj][2]){
					$haplotype3=$pos[$j][$jj][2];$haplotype4=$pos[$j][$jj][3];
				}elsif($AD	eq $pos[$j][$jj][2]){
					$haplotype4=$pos[$j][$jj][2];$haplotype3=$pos[$j][$jj][3];
				}
				if($haplotype1 ne "" and $haplotype2 ne "" and $haplotype3 ne "" and $haplotype4 ne ""){ 
					print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
				}
			}
		}elsif($pos[$j][$jj][4] eq   $pos[$j][$jj][3] and $pos[$j][$jj][5] eq    "H"){
			if ($alt==2 and $heter == 2){
				$haplotype1=$pos[$j][$jj][3];$haplotype2=$pos[$j][$jj][3];
				if($AC	eq $pos[$j][$jj][3]){
					$haplotype3=$pos[$j][$jj][3];$haplotype4=$pos[$j][$jj][2];
				}elsif($AD	eq $pos[$j][$jj][3]){
					$haplotype4=$pos[$j][$jj][3];$haplotype3=$pos[$j][$jj][2];
				}
				if($haplotype1 ne "" and $haplotype2 ne "" and $haplotype3 ne "" and $haplotype4 ne ""){
					print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
				}
			}
		}elsif($pos[$j][$jj][5] eq   $pos[$j][$jj][3] and $pos[$j][$jj][4] eq    "H"){
			if ($alt ==2 and $heter == 2){
				$haplotype3=$pos[$j][$jj][3];$haplotype4=$pos[$j][$jj][3];
				if($AC eq $pos[$j][$jj][3]){
					$haplotype1=$pos[$j][$jj][3];$haplotype2=$pos[$j][$jj][2];
				}elsif($BC eq $pos[$j][$jj][3]){
					$haplotype2=$pos[$j][$jj][3];$haplotype1=$pos[$j][$jj][2];
				}
				if($haplotype1 ne "" and $haplotype2 ne "" and $haplotype3 ne "" and $haplotype4 ne ""){
					print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
				}
			}
		}elsif($pos[$j][$jj][5] eq   $pos[$j][$jj][2] and $pos[$j][$jj][4] eq    "H"){
			if ($ref ==2 and $heter == 2){					
            		$haplotype3=$pos[$j][$jj][2];$haplotype4=$pos[$j][$jj][2];
            		if($AC eq $pos[$j][$jj][2]){
            			$haplotype1=$pos[$j][$jj][2];$haplotype2=$pos[$j][$jj][3];
            		}elsif($BC eq $pos[$j][$jj][2]){
            			$haplotype2=$pos[$j][$jj][2];$haplotype1=$pos[$j][$jj][3];
            		}
					if($haplotype1 ne "" and $haplotype2 ne "" and $haplotype3 ne "" and $haplotype4 ne ""){
						print OUT "$pos[$j][$jj][0]	$pos[$j][$jj][1]	$pos[$j][$jj][2]	$pos[$j][$jj][3]	$pos[$j][$jj][4]	$pos[$j][$jj][5]	$haplotype1	$haplotype2	$haplotype3	$haplotype4\n";
					}
				}
        	}
    }
}
