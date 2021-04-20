#!/bin/bash

#This will run Extract from length 1 to length 5

./Extract < extract.in #First length

for i in {6..16..5}
do
  #Edit and run
  sed -i "s/ $((i-5))/ $i/g" extract.in
  ./Extract < extract.in
done

#Go back to normal
sed -i "s/ 16/ 1/g" extract.in

echo "It is done"
