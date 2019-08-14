# This script parses the first file for the last two BCs.
# It then reads the second file.
# We read the first two BC, we attach to it the whole set of BCs from the first file if they match.

# Input file are count files (or frequency files).

# We need to make sure that the two BCs are unique. They might not be, and if not, we throw them away, and keep track of it, in case we need to correct things.

# However, there are some cases that one of the phasing barcode is not found. 
# Cases:
# BC1-BC2-BC3-BC4
# Phased with:
# BC3-BC4-BC5
# If:
# BC3-BC4, then we use perfect match, and use the highest count BC1-BC2-BC3-BC4 to phase.
# If:
# BC3 is found, but not BC4, then we use the highest count BC1-BC2-BC3 from the previous epoch to phase.
# If:
# BC3 is not found, then search the previous epoch, and phase with the highest BC1-BC2-BC3 from the second-last epoch.
# The highest count is nice, but we want to ensure that the highest count is > 10, since most real previous lineages will have counts over 10. This ensures that we don't match randomly when the counts are close to 1.

open(FILE0,$ARGV[0]) or die; #this is the file with the barcodes from the second-last epoch.

while(my $lines = <FILE0>){
	my @arr = split("\t",$lines);

	my $num_fields = scalar(@arr);

	# How many fields are barcodes?
	my $BC_count = 0;
	my $sum = 0;
	for(my $i = 0;$i < $num_fields;++$i){
		if($arr[$i] =~ m/[ATGCZYXW]+/){
			++$BC_count;
		}
		else{
			$sum += $arr[$i];
		}
	}
	# Ancestry is arr0 to bc_count - 1
	my $ancestry = $arr[0];
	for(my $i = 1;$i < $BC_count-1;++$i){
		$ancestry .= "\t".$arr[$i];
	}

	# Phasing information is arr[bc_count-1]
	my $phase = $arr[$BC_count-1];

	if(!(exists($pre_hash{$phase}))){
		$pre_hash{$phase} = $ancestry;
		$pre_counts{$phase} = $sum;
	}
	else{
		if($sum > $pre_counts{$phase}){
			$pre_hash{$phase} = $ancestry;
			$pre_counts{$phase} = $sum;
		}
	}


}
close(FILE0);

open(FILE1,$ARGV[1]) or die; #this is the file with the barcodes from the last epoch.

while(my $lines = <FILE1>){
	my @arr = split("\t",$lines);

	my $num_fields = scalar(@arr);

	# How many fields are barcodes?
	my $BC_count = 0;
	my $sum = 0;
	for(my $i = 0;$i < $num_fields;++$i){
		if($arr[$i] =~ m/[ATGCZYXW]+/){
			++$BC_count;
		}
		else{
			$sum += $arr[$i];
		}
	}

	# Ancestry is arr0 to bc_count - 2
	my $ancestry = $arr[0];
	for(my $i = 1;$i < $BC_count-2;++$i){
		$ancestry .= "\t".$arr[$i];
	}

	# Phasing information is arr[bc_count-1] and arr[bc_count-2]
	my $phase = $arr[$BC_count-2]."\t".$arr[$BC_count-1];

	if(!(exists($hash{$phase}))){
		$hash{$phase} = $ancestry;
		$counts{$phase} = $sum;
	}
	else{
		if($sum > $counts{$phase}){
			$hash{$phase} = $ancestry;
			$counts{$phase} = $sum;
		}
	}

	# We will not always have complete phasing information (the last barcode may not exist, but we will have matches for the before-last barcode).
	# To solve this, we're going to only use phase information of the before-last barcode, only assuming the most abundant phase.

	$alternate_phases{$arr[$BC_count-2]}{$ancestry} += $sum;
}
close(FILE1);

open(FILE2,$ARGV[2]) or die; # this is the file with counts for the current epoch.

while(my $lines = <FILE2>){
	my @arr = split("\t",$lines);
	if(exists($hash{$arr[0]."\t".$arr[1]})){
		# Case 1
		print $hash{$arr[0]."\t".$arr[1]};
		print "\t";
		print join("\t",@arr);
	}
	elsif(exists($alternate_phases{$arr[0]})){
		# Go through the possible ancestries and use the highest count.
		# Case 2
		my $highest_ancestry;
		my $highest_ancestry_count = 0;
		foreach my $possible_ancestries (keys %{$alternate_phases{$arr[0]}}){
			if($alternate_phases{$arr[0]}{$possible_ancestries} > $highest_ancestry_count){
				$highest_ancestry = $possible_ancestries;
				$highest_ancestry_count = $alternate_phases{$arr[0]}{$possible_ancestries};
			}
		}

		if($highest_ancestry_count >= 10){

			print $highest_ancestry;
			print "\t";
			print join("\t",@arr);
		}
	}
	elsif(exists($pre_hash{$arr[0]})){
		# Case 3
		if($pre_counts{$arr[0]} >= 10){
			print $pre_hash{$arr[0]};
			print "\t";
			print join("\t",@arr);
		}
	}
}