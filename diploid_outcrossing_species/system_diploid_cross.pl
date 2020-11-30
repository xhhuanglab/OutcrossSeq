# perl system_diploid_cross.pl 5 5 100 0.6 0.3 0.35 60 1000 100000 parent1 parent2 demo_diploid demo
# star chromosmoe
# quality threshold we set here as an example,quality below threshold will be filter
# miss rate threshold we set here as an example,missing values greater than the threshold will be filtered
# heterozygosity between 0.3 and (1-0.3)
# start cuttree threshold
# sample number
# min marker number in each window
# window size
# the id of parent1 in vcf file
# the id of parent2 in vcf file
# the input file of the vcf file
# output file's name
$temp1='perl diploid_cross_Genotype_filter.pl '."$ARGV[0] $ARGV[9] $ARGV[10] $ARGV[2] $ARGV[3] $ARGV[11] $ARGV[4]";
print "$temp1\n";
system "$temp1";

$temp2='perl diploid_cross_step1.pl '."$ARGV[8] $ARGV[7] ".'F1-chr'."$ARGV[0]".'.geno '."$ARGV[0] $ARGV[5]";
print "$temp2\n";
system "$temp2";

$temp3='perl diploid_cross_step2.pl '."$ARGV[0] $ARGV[6] $ARGV[5]";
print "$temp3\n";
system "$temp3";

$temp4='perl diploid_cross_step3.pl '."$ARGV[0]";
print "$temp4\n";
system "$temp4";

$temp5='perl diploid_cross_step4.pl '."$ARGV[0] $ARGV[1] demo";
print "$temp5\n";
system "$temp5";
