$qual=$ARGV[1];
$range1=$ARGV[2];
$range2=$ARGV[3];
open OUT1,">ploidy_POS1.vcf";
open OUT2_1,">ploidy_POS2_1.vcf";
open OUT2_2,">ploidy_POS2_2.vcf";
open(IN,"$ARGV[0].vcf")||die "can't open the input file";
chomp(@line=<IN>);
foreach $line(@line){
	if ($line=~/^#CHROM/){
		print OUT1 "$line\n";print OUT2_1 "$line\n";print OUT2_2 "$line\n";
	}elsif($line=~/[GATC]+\s(\d+)/g){
		@lines=split(/\s+/,$line);
		if (length($lines[3]) == 1 and  length($lines[4]) == 1 and $lines[5] >= $qual){
			$line=~s/^\D+//;
			push @array,$line;
			foreach (9..$#lines){
				if($lines[$_]=~/[01]\/[01]\/[01]\/[01]\/[01]\/[01]:(\d+),(\d+):\d+/){
					$alldepth=$1+$2;
					push @alldepth,$alldepth;
				} 
			}
		}
	}
}	
foreach $n(0...1){      
   foreach $i(0...$#alldepth){
      if($i % 2==$n){
         push @{$depth->[$n]},$alldepth[$i];
      }
   }
}
sub quantile{
	my($d,$c,$interval,$quantile1,$quantile2);
	my @depth;
  $d=pop @_; 
  $c=pop @_; 
  @depth=@_;                    
  @depth=sort {$a<=>$b}(@depth); 
  $interval=$#depth;            
  $quantile1=$depth[int($interval*$c)]+($interval*$c-int($interval*$c))*($depth[int($interval*$c)+1]-$depth[int($interval*$c)]);
  $quantile2=$depth[int($interval*$d)]+($interval*$d-int($interval*$d))*($depth[int($interval*$d)+1]-$depth[int($interval*$d)]);
  return($quantile1,$quantile2);
}

foreach $i(0...1){
   ($depth1,$depth2)=&quantile(@{$depth->[$i]},$range1,$range2);
   push @quantile,$depth1;                                      
   push @quantile,$depth2;
}
foreach $line(@array){
	@snp=split(/\s+/,$line);
	foreach $i(9..$#snp){
		if($snp[$i]=~/([01])\/([01])\/([01])\/([01])\/([01])\/([01]):(\d+),(\d+):\d+/){
			$alld=$7+$8;@allele=($1,$2,$3,$4,$5,$6);
			push @alld,$alld;
			$num1=0;$num2=0;
			foreach $j(@allele){
				if ($j==0){
					$num1++;
				}elsif($j==1){
					$num2++;
				}
			}
			$num[$i][0]=$num1;
			$num[$i][1]=$num2;
		}
	}
	foreach $i(0...1){
		if($quantile[2*$i]<$alld[$i] and $alld[$i]<$quantile[2*$i+1]){
			$count++;
		}
	}
	if($count==2){
		if ($num[9][0]==1 and $num[9][1]==5 and $num[10][0]==1 and $num[10][1]==5){#(AT TT TT x AT TT TT)
			print OUT1 "$snp[3]	Chr$line\n";
		}elsif($num[9][0]==5 and $num[9][1]==1 and $num[10][0]==5 and $num[10][1]==1){#(TT TT TA x TT TT TA)
			print OUT1 "$snp[4]	Chr$line\n";
		}elsif($num[9][0]==1 and $num[9][1]==5 and $num[10][1]==6){#(AT TT TT x TT TT TT)
			print OUT2_2 "$snp[3]	Chr$line\n";
		}elsif($num[9][0]==5 and $num[9][1]==1 and $num[10][0]==6){#(AA AA AT x AA AA AA)
			print OUT2_2 "$snp[4]	Chr$line\n";
		}elsif($num[10][0]==5 and $num[10][1]==1 and $num[9][0]==6){#(AA AA AA x AA AA AT)
			print OUT2_1 "$snp[4]	Chr$line\n";
		}elsif($num[10][0]==1 and $num[10][1]==5 and $num[9][1]==6){#(AA AA AA x AA AA AT)
			print OUT2_1 "$snp[3]	Chr$line\n";
		}
	}
	$count=0;
	@alld=();@allele=();
}
close(IN);
