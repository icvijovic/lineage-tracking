# Here we clean counts based on simple rules:
# Must have been observed either:
# 10 times during the whole epoch.
# Or observed in at least two of the time points.

open(FILE,$ARGV[0]) or die;

while(my $lines = <FILE>){
	chomp($lines);
	my @arr = split("\t",$lines);
	my $num_fields = scalar(@arr);
	my $sum = 0;
	my $presence;

	for(my $i = 0;$i < $num_fields;++$i){
		if($arr[$i] !~ m/[ATGCZYXW]+/){
			$sum += $arr[$i];
			if($arr[$i] > 0){
				++$presence;
			}
		}
	}

	#for(my $i = $start;$i < $num_fields;++$i){
	#	$sum += $arr[$i];
	#	if($arr[$i] > 0){
	#		++$presence;
	#	}
	#}

	if($sum >= 10 || $presence >= 2){
		print $lines;
		print "\n";
	}
}