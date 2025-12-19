#!/bin/bash

for i in {4..9}; do ### 1..11
  icx -O3 -o z$i.z plastic11.c -lm

  echo 
done

for i in {4..9}; do 
   echo $i | ./z$i.z & 
   echo 
   sleep 1
   rm z$i.z 
done

