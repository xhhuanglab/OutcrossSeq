$chr = $ARGV[0];
open OUT,">select-$ARGV[1].geno";
$infile="$ARGV[1].geno";
open IN,"$infile" or die "can't open the input file";
while (<IN>){
	chomp();
	$line=$_;
	if ($line=~/Chr(\d+)\s+(\d+)\s+([A T G C])\s+([A T G C])\s+(\d+)\s+(\d+)/){
		if(($6+$5)>1 and $1==$chr and $5/($6+$5)>=0.04 and $5/($6+$5)<=0.28){
			$heter=0;$miss=0;$homo=0;
			@line=split/\s+/,$line;
			foreach (8..$#line){
				if($line[$_] eq "H"){
					$heter++;
				}elsif($line[$_] eq "-"){
					$miss++;
				}else{
					$homo++;
				}
			}
			if (($miss+$homo+$heter)>0 and $miss/($miss+$homo+$heter)<= 0.4 ){
				print OUT "$line\n";
			}
		}
	}elsif($line=~/^#CHROM/){
		print OUT "$line\n";
	}
}
close IN;
close OUT;
