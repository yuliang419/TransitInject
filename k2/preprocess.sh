#! /bin/bash

cd k2mdwarfs

for f in *.txt
do 
	echo $f
	sed -i '' 's/,/ /g' $f
	sed -i '' '1s/^/#/g' $f 
done