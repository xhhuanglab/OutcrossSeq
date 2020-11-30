
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
