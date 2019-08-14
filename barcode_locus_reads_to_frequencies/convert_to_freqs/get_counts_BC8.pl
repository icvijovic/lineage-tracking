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
		print $counts{$bc1}{$bc2}{$bc3}{"7"}{"1"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"7"}{"3"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"7"}{"4"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"7"}{"5"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"7"}{"8"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"7"}{"12"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"8"}{"1"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"8"}{"3"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"8"}{"4"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"8"}{"5"};
		print "\t";
		print $counts{$bc1}{$bc2}{$bc3}{"8"}{"8"};
		print "\n";

		}
	}
}