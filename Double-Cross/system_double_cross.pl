$temp1='perl double_cross_Genotype_filter.pl '."$ARGV[10] $ARGV[6] $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]";
print "$temp1\n";
system "$temp1";

$temp2='perl double_cross_step1.pl '."$ARGV[10]_parent.geno $ARGV[10]_offspring.geno $ARGV[7] $ARGV[4] $ARGV[8] $ARGV[9]";
print "$temp2\n";
system "$temp2";

$temp3='perl double_cross_step2.pl '."$ARGV[10]_parent.geno $ARGV[10]_offspring.geno $ARGV[7] $ARGV[4] $ARGV[8] $ARGV[9]";
print "$temp3\n";
system "$temp3";

$temp4='perl double_cross_step3.pl '."$ARGV[4] $ARGV[5] $ARGV[9] $ARGV[7]";
print "$temp4\n";
system "$temp4";
