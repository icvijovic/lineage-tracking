open(FILE,$ARGV[0]) or die;

while(my $lines = <FILE>){
	my @arr = split("\t",$lines);
	++$counts{$arr[4]}{$arr[5]}{$arr[2]}{$arr[3]};
}

foreach my $bc1 (keys %counts){
	foreach my $bc2 (keys %{$counts{$bc1}}){
		print $bc1;
		print "\t";
		print $bc2;

		print "\t";
		print $counts{$bc1}{$bc2}{"1"}{"1"};
		print "\t";
		print $counts{$bc1}{$bc2}{"1"}{"3"};
		print "\t";
		print $counts{$bc1}{$bc2}{"1"}{"4"};
		print "\t";
		print $counts{$bc1}{$bc2}{"1"}{"5"};
		print "\t";
		print $counts{$bc1}{$bc2}{"1"}{"8"};
		print "\t";
		print $counts{$bc1}{$bc2}{"1"}{"12"};
		print "\t";
		print $counts{$bc1}{$bc2}{"2"}{"1"};
		print "\t";
		print $counts{$bc1}{$bc2}{"2"}{"3"};
		print "\t";
		print $counts{$bc1}{$bc2}{"2"}{"4"};
		print "\t";
		print $counts{$bc1}{$bc2}{"2"}{"5"};
		print "\t";
		print $counts{$bc1}{$bc2}{"2"}{"8"};
		print "\n";
	}
}