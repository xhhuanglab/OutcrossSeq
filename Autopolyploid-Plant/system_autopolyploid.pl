
$temp1='perl  autopolyploid_Genotype_filter_step1.pl '."$ARGV[10] $ARGV[7] $ARGV[2] $ARGV[3]";
print "$temp1\n";
system "$temp1";

$temp2='perl autopolyploid_Genotype_filter_step2.pl '."$ARGV[0] $ARGV[11] $ARGV[8] $ARGV[9]";
print "$temp2\n";
system "$temp2";

$temp3='perl autopolyploid_step1.pl '."$ARGV[0] $ARGV[11]";
print "$temp3\n";
system "$temp3";

$temp4='perl autopolyploid_step2.pl '."select-$ARGV[11] $ARGV[0] $ARGV[4] $ARGV[5] $ARGV[6] $ARGV[7] ";
print "$temp4\n";
system "$temp4";

$temp5='perl autopolyploid_step3.pl '."$ARGV[11]-readsNum.txt R2-select-$ARGV[11] ";
print "$temp5\n";
system "$temp5";

$temp6='perl autopolyploid_step4.pl '." $ARGV[0] R2-select-$ARGV[11]-readsNum.txt R2-select-$ARGV[11].txt ";
print "$temp6\n";
system "$temp6";

$temp7='perl autopolyploid_step5.pl '." $ARGV[0] $ARGV[1] ";
print "$temp7\n";
system "$temp7";