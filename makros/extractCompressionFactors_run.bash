#!/bin/bash
COUNTS_ALL=0
for config in {0,1,2,3,4,5,6,7}; do
#for config in {0,1,3,5}; do
	root -b -q -l extractCompressionFactors_run.C\($config\,5\) > logs/extractCompressionFactors_run_"$config"_5.log &
	root -b -q -l extractCompressionFactors_run.C\($config\,7\) > logs/extractCompressionFactors_run_"$config"_7.log &
	echo "running with generator configuration $config"
	COUNTS_ALL=$((COUNTS_ALL+1))
done
# wait for processes to finish
wait
echo "all ($COUNTS_ALL) jobs done"
