#!/bin/bash

for i in {4..4}; do ### 1..101
 icx -O3 -o z$i.z plastic11.c -lm
#gcc  plastic11.c -g
  echo 
done

for i in {4..4}; do 
   echo $i | ./z$i.z & 
   echo 
   sleep 1
   rm z$i.z 
done

