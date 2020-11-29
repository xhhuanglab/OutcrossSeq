#!/usr/bin/perl-w
$parent1=$ARGV[2];$parent2=$ARGV[3];$parent3=$ARGV[4];$parent4=$ARGV[5];
my ($line,@line,$chr,$ref,$alt,$i,$NumM,$RefDepth,$AlterDepth,$TotalDepth);
open IN,"$ARGV[0].vcf" or die "can't open the input file";
open OUT1,">$ARGV[0]_parent.geno" or die "can't open the output file 1";
open OUT2,">$ARGV[0]_offspring.geno" or die "can't open the output file 2";
while (<IN>){
	chomp();
	$line = $_;
	@line = split/\t+/,$_;
	if($line[0] eq "#CHROM"){
		print OUT1 "$line[0]\t$line[1]\t$line[3]\t$line[4]";
		print OUT2 "$line[0]\t$line[1]\t$line[3]\t$line[4]";
		foreach (9..$#line){
			if ($_ == $parent1-1){
				$parent[0]=$_;
			}elsif($_ == $parent2-1){
				$parent[1]=$_;
			}elsif($_ == $parent3-1){
				$parent[2]=$_;
			}elsif($_ == $parent4-1){
				$parent[3]=$_;
			}else{
				push @offspring,$_;
				print OUT2 "\t$line[$_]";
			}
	    }
		print OUT2 "\n";
		print OUT1 "\t$line[$parent1-1]\t$line[$parent2-1]\t$line[$parent3-1]\t$line[$parent4-1]\n";
    }elsif(length($line[3]) == 1 and length($line[4]) == 1 and $line[5] >= $ARGV[1]){
    	$line[0] =~s/^\D+//;
		$chr = $line[0]+1-1;
    	$ref = $line[3];
    	$alt = $line[4];
		@genotype_parent="";my @genotype;@refrate="";@altrate="";
		$alt_num1=0;$ref_num1=0;$alt_num2=0;$ref_num2=0;
    	foreach $i(@offspring){
    	   $Genotype="-";
		   if ($line[$i]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
			  $RefDepth=$3;$AlterDepth=$4;$TotalDepth=$3+$4;
			  if ($TotalDepth>=1 and $TotalDepth<=500){	
					if ($RefDepth/$TotalDepth>=0.9 ){ 
					   $Genotype=$ref;$ref_num1++;
					}elsif($AlterDepth/$TotalDepth>=0.9){ 
					   $Genotype=$alt;$alt_num1++;
					}elsif ($1==0 and $2==0){ 
					   $Genotype=$ref;$ref_num1++;
					}elsif ($1==1 and $2==1){ 
					   $Genotype=$alt;$alt_num1++;
					}elsif ($AlterDepth/$TotalDepth>=0.4 and  $AlterDepth/$TotalDepth<=0.6){
                       $Genotype="H";
					}elsif ($RefDepth/$TotalDepth>=0.4 and  $RefDepth/$TotalDepth<=0.6){
                       $Genotype="H";
					}elsif ($1==1 and $2==0){
                       $Genotype="H";
					}elsif ($1==0 and $2==1){
                       $Genotype="H"
					}	
				}
			}
			push @genotype,$Genotype;
		}
		foreach $i(@parent){
			$Genotype_parent="-";
			if ($line[$i]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
			  $RefDepth=$3;$AlterDepth=$4;$TotalDepth=$3+$4;
			  if ($TotalDepth>=1 and $TotalDepth<=500){	
					if ($RefDepth/$TotalDepth>=0.9 ){ 
					   $Genotype=$ref;$ref_num2++;
					}elsif($AlterDepth/$TotalDepth>=0.9){ 
					   $Genotype=$alt;$alt_num2++;
					}elsif ($1==0 and $2==0){ 
					   $Genotype=$ref;$ref_num2++;
					}elsif ($1==1 and $2==1){ 
					   $Genotype=$alt;$alt_num2++;
					}elsif ($AlterDepth/$TotalDepth>=0.4 and  $AlterDepth/$TotalDepth<=0.6){
                       $Genotype="H";
					}elsif ($RefDepth/$TotalDepth>=0.4 and  $RefDepth/$TotalDepth<=0.6){
                       $Genotype="H";
					}elsif ($1==1 and $2==0){
                       $Genotype="H";
					}elsif ($1==0 and $2==1){
                       $Genotype="H"
					}	
				}
			}
			$genotype_parent[$i-$parent[0]]=$Genotype;
		}
		if ($ref_num2+$alt_num2==4){
			if ($genotype_parent[0] eq $genotype_parent[1] and $genotype_parent[2] eq $genotype_parent[3] and $genotype_parent[1] eq $genotype_parent[2]){
			}elsif($genotype_parent[0] eq $genotype_parent[1] and $genotype_parent[2] eq $genotype_parent[3] and $genotype_parent[1] ne $genotype_parent[2]){
            }else{
				if ($ref_num1+$alt_num1>0){
					$ref_rate2=$ref_num2/4;$alt_rate2=$alt_num2/4;
					$ref_rate1= sprintf "%0.2f",$ref_num1/($ref_num1+$alt_num1);
					$alt_rate1= sprintf "%0.2f",$alt_num1/($ref_num1+$alt_num1);
					@refrate=($ref_rate1,$ref_rate2);@altrate=($alt_rate1,$alt_rate2);
                    @refrate=(sort {$a <=> $b} @refrate);
                    @altrate=(sort {$a <=> $b} @altrate);
                    $m = $refrate[1]-$refrate[0];$n = $altrate[1]-$altrate[0];
                    if ($m >= 0.2 or $n >= 0.2){
                    }else{
						print OUT1 "$chr\t$line[1]\t$line[3]\t$line[4]";
						print OUT2 "$chr\t$line[1]\t$line[3]\t$line[4]";
						print OUT1 "	$genotype_parent[0]	$genotype_parent[1]	$genotype_parent[2]	$genotype_parent[3]\n";
						foreach $j(0..$#genotype){
							print OUT2 "\t$genotype[$j]";
						}
						print OUT2 "\n";
					}
				}
			}
		}
		
    }	
}
close IN;
close OUT1;
close OUT2;
