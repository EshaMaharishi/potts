#!/bin/bash

echo "numCollagen \t anisotropy"

for i in {1..50}
do 
	dirName="10_23_13/direction${i}"
	./potts $dirName

done;

echo "\n\nfinished all experiments"


