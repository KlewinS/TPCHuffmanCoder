#!/bin/bash
COUNTS_ALL=0
for nCol in {1,2,3,4,5,6,7,8,9,10,11,13,15,17,19}; do
	root -b -q -l generateHuffmanTable_run.C\($nCol\) > logs/generateHuffmanTable_run_"$nCol".log &
	COUNTS_ALL=$((COUNTS_ALL+1))
done
# wait for processes to finish
wait
echo "all ($COUNTS_ALL) jobs done"
