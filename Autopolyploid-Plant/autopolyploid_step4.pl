$chromosome = $ARGV[0];
$file1=$ARGV[1];
$file2= $ARGV[2];
$hang=0;@line = ();%sort=();@array=();@hang_none=();@marker=();@line=();$num=0;
open IN1,"$file1" or die "can't open the input file 1";
open IN2,"$file2" or die "can't open the input file 2";
while (<IN2>){
	$line=$_;
	if ($line=~/^\n/){
		chomp();
		$hang++;
		$sort{$hang}="0";
		push @hang_none,$hang;
	}else{
		chomp();
		$hang++;
		$sort{$hang}="1";
		@line = split/\s+/,$line;
		$marker[$hang][0]=$line[0];$marker[$hang][1]=$line[1];
		$marker[$hang][2]=$line[2];$marker[$hang][3]=$line[3];
		$marker[$hang][4]=$line[4];$marker[$hang][5]=$line[5];
		$marker[$hang][6]=$line[6];$marker[$hang][7]=$line[7];
		foreach $i(8..$#line){
			if($line[$i] ne "H" and $line[$i] ne "-"){
				$geno[$hang]=$line[$i];last;
			}
		}
		foreach $i(8..$#line){
			$marker[$hang][$i*3-16]=$line[$i];
		}
	}
}
unshift @hang_none,0;
$hang=0;
while (<IN1>){
	$array=$_;
	if ($array=~/^\n/){
		chomp();
		$hang++;
		if ($sort{$hang}==0){
		}else{
			print "aa $hang	wrong\n";
		}
	}else{
		chomp();
		$hang++;
		@array=split/\s+/,$array; 
		if ($sort{$hang}==1 and $marker[$hang][0] eq $array[0] and $marker[$hang][1] eq $array[1]){
			foreach $i(8..$#array){
				@read=split/\,/,$array[$i];
				$marker[$hang][3*$i-15]=$read[0];
				$marker[$hang][3*$i-14]=$read[1];
			}
		}else{
			print "bb $hang	wrong\n";
		}
	}
}
$out="Chr$chromosome".'-Imputation_marker.txt';
open OUT,">$out" or die "can't open the output file";
foreach $i(1..$hang){
	if ($sort{$i}==0 and $i >1){
			$num++;
		foreach $jj(8..$#array){
			$miss = 0;$Heterozygote=0;$reads_num=0;$reads_num1=0;$reads_num2=0;
			foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
				if ($marker[$ii][3*$jj-16] eq "-"){
					$miss++;
				}elsif($marker[$ii][3*$jj-16] eq "H"){
					$Heterozygote++;
				}
				if ($marker[$ii][3*$jj-16] ne "-" and $marker[$ii][3*$jj-16] ne "H"){
					$reads_num=$reads_num+$marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14];
					$reads_num1=$reads_num1+$marker[$ii][3*$jj-15];
					$reads_num2=$reads_num2+$marker[$ii][3*$jj-14];
				}
				if($marker[$ii][3*$jj-16] eq "H"){
					$reads_num1=$reads_num1+$marker[$ii][3*$jj-15];
				}
			}
			$marker_num = $hang_none[$num]-1-($hang_none[$num-1]+1)+1;
			if ($marker_num>=10 and $marker_num <=20){
				if($marker_num-$miss >= 8){
					if($Heterozygote >= 3){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= "H";
						}
					}elsif ($Heterozygote >= 1 and $Heterozygote <= 2){
						if($reads_num > 0 and $reads_num1/$reads_num <= 0.03){
							foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
								if ($marker[$ii][3*$jj-16] eq "H"){
									$marker[$ii][3*$jj-16]= "H";
								}else{
									$marker[$ii][3*$jj-16]= $geno[$ii];
								}
							}
						}elsif($reads_num > 0 and $reads_num1/$reads_num <= 0.05 and $reads_num1/$reads_num > 0.03){
							foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
								if ($marker[$ii][3*$jj-16] eq "H"){
									$marker[$ii][3*$jj-16]= "H";
								}elsif(($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14]) <= 6){
									$marker[$ii][3*$jj-16]= "H";
								}else{
									$marker[$ii][3*$jj-16]= $geno[$ii];
								}
							}
								
						}else{
							foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
								if (($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14]) < 10){
									$marker[$ii][3*$jj-16]= "H";
								}
							}
						}
					}elsif (($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.95 and $reads_num/($marker_num-$miss-$Heterozygote)>=7){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}	
					}elsif($reads_num1	==0 and $reads_num/($marker_num-$miss-$Heterozygote)>=6){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif(($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.8 ){
						if( $reads_num1	==0 and $reads_num/($marker_num-$miss-$Heterozygote)<=4){
							foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
								if($marker[$ii][3*$jj-14] >=5){
									$marker[$ii][3*$jj-16]= $geno[$ii];
								}else{
									$marker[$ii][3*$jj-16]= "H";
								}
							}
						}else{
							foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
								if ($marker[$ii][3*$jj-16] eq "H"){
									$marker[$ii][3*$jj-16]= "H";
								}else{
									$marker[$ii][3*$jj-16]= $geno[$ii];
								}
							}
						}
					}else{
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if ($marker[$ii][3*$jj-16] eq "H"){
								$marker[$ii][3*$jj-16]= "H";
							}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 5){
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}else{
								$marker[$ii][3*$jj-16]= "H";
							}
						}
					}
				}else{
					if ($Heterozygote >= 3){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= "H";
						}
					}elsif($Heterozygote >= 1 and $Heterozygote < 3){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if ($marker[$ii][3*$jj-16] eq "H"){
								$marker[$ii][3*$jj-16]= "H";
							}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 7){
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}else{
								$marker[$ii][3*$jj-16]= "H";
							}
						}
					}elsif(	$reads_num >= 20 and $reads_num1==0 ){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}else{
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 6){
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}else{
								$marker[$ii][3*$jj-16]= "H";
							}
						}
					}
				}
			}elsif($marker_num>=21 and $marker_num<=50){
				if ($miss/$marker_num < 0.6){
					if ($Heterozygote/($marker_num-$miss) >= 0.4 or $Heterozygote >= 5){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= "H";
						}
					}elsif($Heterozygote/($marker_num-$miss) >= 0.2){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if ($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 15 and $marker[$ii][3*$jj-16] ne "H"){
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}else{
								$marker[$ii][3*$jj-16]= "H";
							}
						}
					}elsif (($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.95 and $reads_num/($marker_num-$miss-$Heterozygote)>=7){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif($reads_num1	==0 and $reads_num/($marker_num-$miss-$Heterozygote)>=5){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif(($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.8 ){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if ($marker[$ii][3*$jj-16] eq "H"){
								$marker[$ii][3*$jj-16]= "H";
							}else{
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}
						}
					}
				}elsif($Heterozygote>0 and $Heterozygote/($marker_num-$miss) >= 0.4 or $Heterozygote >= 5){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						$marker[$ii][3*$jj-16]= "H";
					}
				}elsif($Heterozygote >= 1 and $Heterozygote < 5){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						if ($marker[$ii][3*$jj-16] eq "H"){
							$marker[$ii][3*$jj-16]= "H";
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 8){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
							$marker[$ii][3*$jj-16]= "H";
						}else{
							$marker[$ii][3*$jj-16]= "H";
						}
					}
				}elsif(	$reads_num >= 50 and $reads_num1==0 ){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						$marker[$ii][3*$jj-16]= $geno[$ii];
					}
				}else{
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						if($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 6){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
							$marker[$ii][3*$jj-16]= "-";
						}else{
							$marker[$ii][3*$jj-16]= "H";
						}
					}
				}
			}elsif($marker_num>=51 and $marker_num<=80){
				if ($miss/$marker_num < 0.7){
					if ($Heterozygote/($marker_num-$miss) >= 0.4 or $Heterozygote >= 8){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= "H";
						}
					}elsif($Heterozygote/($marker_num-$miss) >= 0.2){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if ($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 15 and $marker[$ii][3*$jj-16] ne "H"){
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}else{
								$marker[$ii][3*$jj-16]= "H";
							}
						}
					}elsif (($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.95 and $reads_num/($marker_num-$miss-$Heterozygote)>=7){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif($reads_num1	==0 and $reads_num/($marker_num-$miss-$Heterozygote)>=5){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif(($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.8 ){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
								$marker[$ii][3*$jj-16]= "-";
							}
						}
					}
				}elsif($Heterozygote>0 and $Heterozygote/($marker_num-$miss) >= 0.4 or $Heterozygote >= 5){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						$marker[$ii][3*$jj-16]= "H";
					}
				}elsif($Heterozygote >= 1 and $Heterozygote < 5){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						if ($marker[$ii][3*$jj-16] eq "H"){
							$marker[$ii][3*$jj-16]= "H";
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 8){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
							$marker[$ii][3*$jj-16]= "-";
						}else{
							$marker[$ii][3*$jj-16]= "H";
						}
					}
				}elsif(	$reads_num >= 80 and $reads_num1==0 ){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						$marker[$ii][3*$jj-16]= $geno[$ii];
					}
				}else{
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						if($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 6){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
							$marker[$ii][3*$jj-16]= "-";
						}else{
							$marker[$ii][3*$jj-16]= "H";
						}
					}
				}
			}elsif($marker_num>=81){
				if ($miss/$marker_num < 0.8){
					if ($Heterozygote/($marker_num-$miss) >= 0.4 or $Heterozygote >= 10){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= "H";
						}
					}elsif($Heterozygote/($marker_num-$miss) >= 0.2){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if ($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 15 and $marker[$ii][3*$jj-16] ne "H"){
								$marker[$ii][3*$jj-16]= $geno[$ii];
							}else{
								$marker[$ii][3*$jj-16]= "H";
							}
						}
					}elsif (($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.95 and $reads_num/($marker_num-$miss-$Heterozygote)>=7){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif($reads_num1	==0 and $reads_num/($marker_num-$miss-$Heterozygote)>=5){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}
					}elsif(($marker_num-$miss-$Heterozygote)/($marker_num-$miss) >= 0.8 ){
						foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
							if($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
								$marker[$ii][3*$jj-16]= "-";
							}
						}
					}
				}elsif($Heterozygote>0 and $Heterozygote/($marker_num-$miss) >= 0.4 or $Heterozygote >=7){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						$marker[$ii][3*$jj-16]= "H";
					}
				}elsif($Heterozygote >= 1 and $Heterozygote < 7){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						if ($marker[$ii][3*$jj-16] eq "H"){
							$marker[$ii][3*$jj-16]= "H";
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 8){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
							$marker[$ii][3*$jj-16]= "-";
						}else{
							$marker[$ii][3*$jj-16]= "H";
						}
					}
				}elsif(	$reads_num >= 100 and $reads_num1==0 ){
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						$marker[$ii][3*$jj-16]= $geno[$ii];
					}
				}else{
					foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
						if($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] >= 6){
							$marker[$ii][3*$jj-16]= $geno[$ii];
						}elsif($marker[$ii][3*$jj-15]+$marker[$ii][3*$jj-14] < 3){
							$marker[$ii][3*$jj-16]= "-";
						}else{
							$marker[$ii][3*$jj-16]= "H";
						}
					}
				}
			}
		}
		foreach $ii($hang_none[$num-1]+1..$hang_none[$num]-1){
			print OUT "$marker[$ii][0]	$marker[$ii][1]	$marker[$ii][2]	$marker[$ii][3]";
			foreach $jj(8..$#array){
				print OUT "\t$marker[$ii][3*$jj-16]";
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
}
close IN1;
close IN2;
close OUT;	
$outfile="Chr$chromosome".'-Imputation_marker_sort.txt';
open IN3,"$out" or die "can't open the file $out";
open IN4,"sample.list" or die "can't open the file sample.list";
open OUT2,">$outfile" or die "can't open the out file ";
print OUT2 "#CHROM	Position	Ref	Alt";
while (<IN4>){
	chomp();
	print OUT2 "\t$_";
}
print OUT2 "\n";
our @pos;
while (<IN3>){
	chomp();
	$snp=$_;
	if ($snp=~/^Chr(\d+)\s+(\d+)/){
		if($chromosome==$1){
			$pos=$2;
			$snp{$pos}=$snp;
			push @pos,$pos;
		}
	}
}
@sortted_pos = sort {$a <=> $b}@pos;
foreach (@sortted_pos){
	print OUT2 "$snp{$_}\n";
}
close OUT2 ;
close IN4;
close IN3;