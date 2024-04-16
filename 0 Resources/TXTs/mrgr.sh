#!/bin/bash

counter=1
separator=$(printf '%.0s-' {1..64})
for file in *.txt
do
  echo "$separator" >> merge.txt
  echo "${counter}CHCK" >> merge.txt
  echo "$separator" >> merge.txt
  cat "$file" >> merge.txt
  ((counter++))
done