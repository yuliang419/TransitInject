#! /bin/bash

cd k2mdwarfs
sed -i '' 's/,/ /g' *.txt
for f in *.txt; do sed -i '' '1s/^/#/g' $f; done