open(FILE,$ARGV[0]) or die;

while(my $lines = <FILE>){
	my @arr = split("\t",$lines);
	++$counts{$arr[4]}{$arr[5]}{$arr[6]}{$arr[2]}{$arr[3]};
}

foreach my $bc1 (keys %counts){
	foreach my $bc2 (keys %{$counts{$bc1}}){
		foreach my $bc3 (keys %{$counts{$bc1}{$bc2}}){
			print $bc1;
			print "\t";
			print $bc2;
			print "\t";
			print $bc3;

			print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"5"}{"1"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"5"}{"3"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"5"}{"4"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"5"}{"5"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"5"}{"7"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"5"}{"12"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"6"}{"1"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"6"}{"3"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"6"}{"4"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"6"}{"5"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"6"}{"7"};
		print "\n";

		}
	}
}