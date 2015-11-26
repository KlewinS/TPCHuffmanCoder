#!/bin/bash
COUNTS_ALL=0
for config in {0,1,2,3,4,5,6,7}; do
	for nCol in {5,7,9}; do
		root -b -q -l generateHuffmanTable_run.C\($nCol\,$config\) > logs/generateHuffmanTable_run_"$nCol"_"$config".log &
		COUNTS_ALL=$((COUNTS_ALL+1))
	done
done
# wait for processes to finish
wait
echo "all ($COUNTS_ALL) jobs done"
