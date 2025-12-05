#!/bin/sh

export OMP_NUM_THREADS=8
DirName=GMSfiles
for name in `ls -1 $DirName`
do
	echo "$name"
	./sparsePOP $DirName/$name
done
