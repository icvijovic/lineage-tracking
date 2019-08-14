# We open the file that will map the frequencies.

open(MAP,$ARGV[0]) or die;
my @readcounts;

my $n = 0;
my $f = 5;

my $number_of_files = 0;
my $skip_file = 0;
my $max_skip = $f-2;

while(my $lines = <MAP>){
	chomp($lines);
	if(-f $lines){
		
		open(FILE,$lines) or die;
		if($skip_file < $max_skip){
			++$skip_file;
			close(FILE);
			next;
		}

		while(my $lines2 = <FILE>){
			chomp($lines2);
			my @arr = split("\t",$lines2);
			my $individual = $arr[0]."\t".$arr[1]."\t".$arr[2]."\t".$arr[3]."\t".$arr[4];

			$pops{$individual}[$n] += $arr[$f];
			$pops{$individual}[$n+1] += $arr[$f + 1];
			$pops{$individual}[$n+2] += $arr[$f + 2];
			$pops{$individual}[$n+3] += $arr[$f + 3];
			$pops{$individual}[$n+4] += $arr[$f + 4];
			$pops{$individual}[$n+5] += $arr[$f + 5];
			$pops{$individual}[$n+6] += $arr[$f + 6];
			$pops{$individual}[$n+7] += $arr[$f + 7];
			$pops{$individual}[$n+8] += $arr[$f + 8];
			$pops{$individual}[$n+9] += $arr[$f + 9];
			$pops{$individual}[$n+10] += $arr[$f + 10];

			$readcounts[$n] += $arr[$f];
			$readcounts[$n+1] += $arr[$f + 1];
			$readcounts[$n+2] += $arr[$f + 2];
			$readcounts[$n+3] += $arr[$f + 3];
			$readcounts[$n+4] += $arr[$f + 4];
			$readcounts[$n+5] += $arr[$f + 5];
			$readcounts[$n+6] += $arr[$f + 6];
			$readcounts[$n+7] += $arr[$f + 7];
			$readcounts[$n+8] += $arr[$f + 8];
			$readcounts[$n+9] += $arr[$f + 9];
			$readcounts[$n+10] += $arr[$f + 10];
		}

		close(FILE);
		$n += 11;
		$f += 1;

		++$number_of_files;
	}
}

my $offset = 330;
print "POP	BC1	BC2	BC3	BC4";

for(my $i = 0;$i < $number_of_files;++$i){
	print "\t";
	print $offset;
	print "\t";
	print ($offset+10);
	print "\t";
	print ($offset+20);
	print "\t";
	print ($offset+30);
	print "\t";
	print ($offset+40);
	print "\t";
	print ($offset+50);
	print "\t";
	print ($offset+60);
	print "\t";
	print ($offset+70);
	print "\t";
	print ($offset+80);
	print "\t";
	print ($offset+90);
	print "\t";
	print ($offset+100);

	$offset += 110;
}
print "\n";

foreach my $individual (keys %pops){
	print $individual;
	for(my $i = 0;$i < scalar(@readcounts);++$i){
		print "\t";
		print ($pops{$individual}[$i]/$readcounts[$i]);
			
	}
	print "\n";
}
for(my $i = 0;$i < scalar(@readcounts);++$i){
	print STDERR $i;
	print STDERR "\t";
	print STDERR $readcounts[$i];
	print STDERR "\n";		
}