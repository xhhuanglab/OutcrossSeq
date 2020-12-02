use List::Util qw/max min/;
open IN1,"ploidy_POS1.vcf" or die "can't open the input file";
open IN2,"ploidy_POS2_1.vcf" or die "can't open the input file";
open IN3,"ploidy_POS2_2.vcf" or die "can't open the input file";

while (<IN1>){
	chomp();
	if($_=~/([A T G C])	Chr(\d+)\s+(\d+)/){
		if ($ARGV[0] == $2){
			push @key1,"$2	$3";
			$a="$2	$3";
			$less_geno{$a}=$1;
		}elsif($ARGV[0] < $1){
			last;
		}
	}
}
while (<IN2>){
	chomp();
	if($_=~/([A T G C])	Chr(\d+)\s+(\d+)/){
		if ($ARGV[0] == $2){
			push @key2,"$2	$3";
			$a="$2	$3";
			$less_geno{$a}=$1;
		}elsif($ARGV[0] < $1){
			last;
		}
	}
}
while (<IN3>){
	chomp();
	if($_=~/([A T G C])	Chr(\d+)\s+(\d+)/){
		if ($ARGV[0] == $2){
			push @key3,"$2	$3";
			$a="$2	$3";
			$less_geno{$a}=$1;
		}elsif($ARGV[0] < $1){
			last;
		}
	}
}
open IN4,"$ARGV[1].vcf" or die "can't open the input file 4";
open OUT1,">$ARGV[1].geno";
open OUT2,">$ARGV[1]-readsNum.txt";
while (<IN4>){
	chomp();
	$line = $_;
	if ($line=~/^#CHROM/){
		@line = split/\s+/,$line;
		print OUT1 "$line[0]\t$line[1]\t$line[3]\t$line[4]\tLessDepth\tMoreDepth\tHeterozygousRate\tHomozygous-Rate";
		print OUT2 "$line[0]\t$line[1]\t$line[3]\t$line[4]\tLessDepth\tMoreDepth\tHeterozygousRate\tHomozygous-Rate";
		foreach $i(9..$#line){
			if ($line[$i] ne "$ARGV[2]" and $line[$i] ne "$ARGV[3]"){
				push @num,$i;
				print OUT1 "\t$line[$i]";
				print OUT2 "\t$line[$i]";
			}
		}
		print OUT1 "\n";print OUT2 "\n";
	}elsif($line=~/[GATC]+\s(\d+)/g){
		@line = split/\s+/,$line;
		$line[0] =~s/^\D+//;$chr=$line[0];
		if ($chr ==$ARGV[0]){
			@line = split/\s+/,$line;$pos=$line[1];$key = "$chr	$pos";
			if ( grep { $_ eq $key } @key1 or grep { $_ eq $key } @key2 or grep { $_ eq $key } @key3) {
				if ($line[5] >= 100 and $line[7] ne "LowQual" and length($line[3])==1 and length($line[4])==1){
					$LessDepth=0;$MoreDepth=0;$miss_num=0;$refNum=0;$heter=0;$altNum=0;$ref = $line[3];$alt=$line[4];
					@Genotype=();@less_geno_num=();@more_geno_num=();
					foreach $i(@num){
						$Genotype="-";$less_geno_num=0;$more_geno_num=0;
						if ($line[$i]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
							if ($less_geno{$key} eq $ref){
								$less_geno_num=$3;$more_geno_num=$4;$LessDepth=$3+$LessDepth;$MoreDepth=$4+$MoreDepth;
							}elsif($less_geno{$key} eq $alt){
								$less_geno_num=$4;$more_geno_num=$3;$LessDepth=$4+$LessDepth;$MoreDepth=$3+$MoreDepth;
							}
							if(($3+$4)>=3 and ($3+$4)<=50){
								if ($3/($3+$4)>=0.98 ){ 
									$Genotype=$ref;
								}elsif($4/($3+$4)>=0.98){ 
									$Genotype=$alt;
								}elsif ($4/($3+$4)>=0.04 and  $4/($3+$4)<=0.28){
									$Genotype="H";
								}elsif ($3/($3+$4)>=0.04 and  $3/($3+$4)<=0.28){
									$Genotype="H";
								}elsif($1!=$2){
									$Genotype="H";
								}elsif($1==$2 and $1==0){
									$Genotype=$ref;
								}elsif($1==$2 and $1==1){
									$Genotype=$alt;
								}
							}
							if (($3+$4)<=50){
								if ($less_geno_num >= 1){
									$Genotype="H";
								}
							}
						}
						if ($Genotype eq "-"){
							$miss_num++;
						}elsif($Genotype eq $ref){
							$refNum++;
						}elsif($Genotype eq $alt){
							$altNum++;
						}elsif($Genotype eq "H"){
							$heter++;
						}
						push @Genotype,$Genotype;
						push @less_geno_num,$less_geno_num;
						push @more_geno_num,$more_geno_num;
					}
					if($#line-8-$miss_num > 0){
						$rate_ref=sprintf "%0.4f",($refNum/($#line-8-$miss_num));
						$rate_alt=sprintf "%0.4f",($altNum/($#line-8-$miss_num));
						$rate_heter=sprintf "%0.4f",($heter/($#line-8-$miss_num));
					}else{
						$rate_ref=0;$rate_alt=0;$rate_heter=0;
					}
					if($less_geno{$key} eq $ref){
						$heter_rate=$rate_heter;
						$homo_rate=$rate_alt;
					}else{
						$heter_rate=$rate_heter;
						$homo_rate=$rate_ref;
					}
					print OUT1 "Chr$chr\t$line[1]\t$line[3]\t$line[4]\t$LessDepth\t$MoreDepth\t$heter_rate\t$homo_rate";
					print OUT2 "Chr$chr\t$line[1]\t$line[3]\t$line[4]\t$LessDepth\t$MoreDepth\t$heter_rate\t$homo_rate";
					foreach (0..$#Genotype){
						print OUT1 "\t$Genotype[$_]";
						print OUT2 "\t$less_geno_num[$_],$more_geno_num[$_]";
					}
					print OUT1 "\n";
					print OUT2 "\n";
				}
			}
		}
	}
}
close IN1;
close IN2;
close IN3;
close IN4;
close OUT2;
close OUT1;

			