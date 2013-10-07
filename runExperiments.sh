#!/bin/bash

echo "numCollagen \t anisotropy"

for i in {1..50}
do 
	dirName="output${i}"
	./potts $dirName $i 10 
done;

echo "\n\nfinished all experiments"


