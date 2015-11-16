#!/bin/bash
COUNTS_ALL=0
for bits in {10,12}; do
	root -b -q -l extractProbabilities_run.C\($bits\) > logs/extractProbabilities_run_"$bits".log &
	COUNTS_ALL=$((COUNTS_ALL+1))
done
# wait for processes to finish
wait
echo "all ($COUNTS_ALL) jobs done"
