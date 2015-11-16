#!/bin/bash
COUNTS_ALL=0
for rate in {1,3,5,7,9,11,13,15,17,19}; do
	root -b -q -l extractCompressionFactors_fromData_run.C\($rate\) > logs/extractCompressionFactors_fromData_run_"$rate".log &
	echo "running with $rate collisions table"
	COUNTS_ALL=$((COUNTS_ALL+1))
done
# wait for processes to finish
wait
echo "all ($COUNTS_ALL) jobs done"
