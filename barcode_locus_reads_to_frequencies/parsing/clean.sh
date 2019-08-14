clean () (
	gunzip -c $* |
	awk -F : 'BEGIN {OFS="\t"}
		NR%4==1 {idx=$10}
		NR%4==2 {print idx, $0}
	' |
	tr '+' '\t'
)

paste <(clean $1) <(clean $2 | awk '{print $3}') | grep -v -P 'G{10}(	|$)'