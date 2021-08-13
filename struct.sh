#!/bin/bash


for (( a = 2; a < 32769; a*=2 ))
do
./struct.out "$a"
done