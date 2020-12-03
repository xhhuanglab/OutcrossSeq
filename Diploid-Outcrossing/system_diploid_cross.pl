
$temp1='perl diploid_cross_Genotype_filter.pl '."$ARGV[0] $ARGV[8] $ARGV[9] $ARGV[1] $ARGV[2] $ARGV[10] $ARGV[3]";
print "$temp1\n";
system "$temp1";

$temp2='perl diploid_cross_step1.pl '."$ARGV[7] $ARGV[6] ".'F1-chr'."$ARGV[0]".'.geno '."$ARGV[0] $ARGV[4]";
print "$temp2\n";
system "$temp2";

$temp3='perl diploid_cross_step2.pl '."$ARGV[0] $ARGV[5] $ARGV[4]";
print "$temp3\n";
system "$temp3";

$temp4='perl diploid_cross_step3.pl '."$ARGV[0]";
print "$temp4\n";
system "$temp4";
