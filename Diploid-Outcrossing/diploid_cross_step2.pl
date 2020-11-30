#perl diploid_cross_step2.pl 1 184 0.35
$chr = $ARGV[0];
$command='ls Chr'."$chr".'_*-'."$ARGV[2]".'.kinf > Chr'."$chr".'-kinf.list';
print "$command\n";
system ("$command");
open IN1,"Chr$chr-kinf.list" or die "can't open the input file1";
while (<IN1>){
	chomp();
	if ($_=~/Chr(\d+)\_(\d+)_(\d+)/){
		$filename{$2}=$_;
		$pos{$2}=$3;
		push @file_num,$2;
	}
}
@sort_file_num=sort { $a <=> $b } @file_num;
$sample_num=$ARGV[1];
$min_sample = $sample_num/4*0.7;
$max_sample = $sample_num/4*1.3; 
$out_geno="Chr$chr".'.genotype';
open OUT_geno,">$out_geno" or die "can't open the OUT_geno file";
foreach $sort_file_num(@sort_file_num){
	$group_num = $sample_num;$cuttree = $ARGV[2];@type_num_ture_sample_num=();$type_num_ture=0;
	while ($group_num >= $max_sample){
		$input_kin = 'Chr'."$chr".'_'."$sort_file_num".'_'."$pos{$sort_file_num}";
		$temp='Rscript diploid_cross_cluster.R '."$input_kin"." $cuttree";
		print "$temp\n";
		system ("$temp");
		@species=();@kin_name=();@sample_name1=();$group_num = 0;@type_num=();@species=();@sample_name=();
		$input_genotype='Chr'."$chr".'_'."$sort_file_num".'_'."$pos{$sort_file_num}".'_'."$cuttree".'.genotype';
		open IN ,"$input_genotype" or die "can't open the genotype $chromosome $sort_file_num";
		while (<IN>){
			chomp();
			@line = split/\s+/,$_;
			$type{$line[0]}=$line[1];
			push @sample_name,$line[0];
			if ( grep { $_== $line[1] } @species){
			}else{
				push @species,$line[1];
			}
		}
		@sort_species=sort { $b <=> $a } @species;
		foreach $sample_name(@sample_name){
			foreach $i(1..$sort_species[0]){
				if ($type{$sample_name} == $i){
					$type_num[$i]++;
				}
			}
		}
		$cuttree=sprintf "%0.2f",$cuttree-0.01;
		foreach $i(1..$sort_species[0]){
			if ($type_num[$i] >= $min_sample and $type_num[$i]<=$max_sample){
				$type_num_ture++;$type_num_ture_sample_num=0;
				push @type_num_ture_sample_num,$type_num[$i];
				foreach $j(0..$#sample_name){
					if($type{$sample_name[$j]} == $i){
						$type_num_ture_sample_num++;
						$sample_type[$type_num_ture][$type_num_ture_sample_num]="$sample_name[$j]:$type_num_ture";
					}
				}
			}elsif($type_num[$i] > $max_sample){
				$group_num = $group_num+$type_num[$i];
				foreach $j(0..$#sample_name){
					if ($type{$sample_name[$j]} eq $i){
						push @sample_name1,$sample_name[$j];
					}
				}
			}
		}
		if ($group_num > 0 ){
			$outfile='Chr'."$chr".'_'."$sort_file_num".'_'."$pos{$sort_file_num}".'-'."$cuttree".'.kinf';
			open OUT,">$outfile";
			$hang=0;$match=0;
			open IN1,"$input_kin-$ARGV[2].kinf" or die "can't open the in file";	
			while (<IN1>){
				chomp();
				$hang++;
				if ($hang==1){
					@kin = split/\s+/,$_;
					foreach $m(0..$#kin){
						foreach $n(@sample_name1){
							if ($n eq $kin[$m]){
								$match++;
								if ($match==1){
									push @kin_name,$m;print OUT "$n";
								}else{
									push @kin_name,$m;print OUT "\t$n";
								}
							}
						}
					}
					print OUT "\n";
				}else{
					@kin = split/\s+/,$_;
					if ( grep { $_==$hang-2 } @kin_name){
						foreach $m(0..$#kin_name){
							if($m==0){
								print OUT "$kin[$kin_name[$m]]";
							}else{
								print OUT "\t$kin[$kin_name[$m]]";
							}
						}print OUT "\n";
					}
				}
			}
		}else{
		$group_num=0;
		}
	}
	print OUT_geno "$chr	$sort_file_num	$pos{$sort_file_num}";
	foreach $i(1..$type_num_ture){
		foreach $j(1..$type_num_ture_sample_num[$i-1]){
			print OUT_geno "	$sample_type[$i][$j]";
		}	
	}
	print OUT_geno "\n";
}	
