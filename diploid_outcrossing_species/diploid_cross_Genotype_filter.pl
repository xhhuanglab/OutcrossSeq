#!/usr/bin/perl-w
$qual = $ARGV[3];$chr = $ARGV[0];$miss_rate = $ARGV[4];
open IN,"$ARGV[5].vcf" or die "can't open the input file";
open OUT2,">F1-chr$ARGV[0].geno" or die "can't open the output file 1";
open OUT1,">Parents-chr$ARGV[0].geno" or die "can't open the output file 2";
while (<IN>){
	chomp();
	$line = $_;
	@line = split/\t+/,$_;
	if($line[0] eq "#CHROM"){
		print OUT1 "$line[0]\t$line[1]\t$line[3]\t$line[4]";
		print OUT2 "$line[0]\t$line[1]\t$line[3]\t$line[4]";
		foreach (9..$#line){
			if($line[$_] ne $ARGV[1] and $line[$_] ne $ARGV[2]){
				print OUT2 "\t$line[$_]";
			}else{
				push @parent,$_;
				print OUT1 "\t$line[$_]";
			}
	    }
	    print OUT1 "\n";
		print OUT2 "\n";
	}elsif($line=~/[GATC]+\s(\d+)/g ){
		$line[0] =~s/^\D+//;
		@genotype=();$miss=0;
		if($line[0]==$chr and length($line[3]) == 1 and length($line[4]) == 1 and $line[5] >= $qual){
			$ref = $line[3];$alt = $line[4];
			foreach $i(0..$#parent){
				$parent_geno[$i]="-";
				if ($line[$parent[$i]]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
					$RefDepth=$3;$AlterDepth=$4;$TotalDepth=$3+$4;
					if ($TotalDepth>=10 and $TotalDepth<=1000){	
						if ($RefDepth/$TotalDepth>=0.9 ){ 
							$parent_geno[$i]=$ref;
						}elsif($AlterDepth/$TotalDepth>=0.9){ 
							$parent_geno[$i]=$alt;
						}elsif ($AlterDepth/$TotalDepth>=0.4 and  $AlterDepth/$TotalDepth<=0.6){
							$parent_geno[$i]="H";
						}elsif ($RefDepth/$TotalDepth>=0.4 and  $RefDepth/$TotalDepth<=0.6){
							$parent_geno[$i]="H";
						}elsif ($1==0 and $2==0){ 
							$parent_geno[$i]=$ref;
						}elsif ($1==1 and $2==1){ 
							$parent_geno[$i]=$alt;$NumM+=1;
						}elsif ($1==1 and $2==0){
							$parent_geno[$i]="H";
						}elsif ($1==0 and $2==1){
							$parent_geno[$i]="H";
						}
					}
				}
			}
			
			if ($parent_geno[0] eq "H" or $parent_geno[1] eq "H"){
				@genotype=();$miss=0;$ref = $line[3]; $alt = $line[4];$heter=0;$homo=0;
				foreach $i(9..$#line){
					if ($i !=  $parent[0] and $i != $parent[1]){
						$Genotype="-";
						if ($line[$i]=~/(\d)\/(\d)\:(\d+)\,(\d+)\:/){
							$RefDepth=$3;$AlterDepth=$4;$TotalDepth=$3+$4;
							if ($TotalDepth>=1 and $TotalDepth<=500){
								if ($RefDepth/$TotalDepth>=0.9 ){ 
									$Genotype=$ref;
								}elsif($AlterDepth/$TotalDepth>=0.9){ 
									$Genotype=$alt;
								}elsif ($AlterDepth/$TotalDepth>=0.4 and  $AlterDepth/$TotalDepth<=0.6){
									$Genotype="H";
								}elsif ($RefDepth/$TotalDepth>=0.4 and  $RefDepth/$TotalDepth<=0.6){
									$Genotype="H";
								}elsif ($1==0 and $2==0){ 
									$Genotype=$ref;
								}elsif ($1==1 and $2==1){ 
									$Genotype=$alt;$NumM+=1;
								}elsif ($1==1 and $2==0){
									$Genotype="H";
								}elsif ($1==0 and $2==1){
									$Genotype="H"
								}
							}
						}
						push @genotype,$Genotype;
						if ($Genotype eq "-"){
							$miss++;
						}elsif($Genotype eq "H"){
							$heter++;
						}elsif($Genotype eq $ref or $Genotype eq $alt){
							$homo++;
						}
					}
				}
				if ($miss/($#genotype+1)<=$miss_rate){
					if ($heter/($homo+$heter)>=$ARGV[6] and $heter/($homo+$heter)<=1-$ARGV[6]){
						print OUT1 "Chr$line[0]\t$line[1]\t$line[3]\t$line[4]\t$parent_geno[0]\t$parent_geno[1]\n";
						print OUT2 "Chr$line[0]\t$line[1]\t$line[3]\t$line[4]";
						foreach (@genotype){
							print OUT2 "\t$_";
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
					
							
				
		
