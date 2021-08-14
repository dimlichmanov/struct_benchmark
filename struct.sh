#!/bin/bash


for (( a = 1; a < 262145; a*=2 ))
do
./struct.out "$a"
done