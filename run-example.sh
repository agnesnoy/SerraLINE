#!/bin/bash

#This will run Extract from length 1 to length 5

./Extract < extract.in #First length

for i in {2..5}
do
  #Edit and run
  sed -i "s/$((i-1))/$i/g" extract.in
  ./Extract < extract.in
done

echo "It is done"
