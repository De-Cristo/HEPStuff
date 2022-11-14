#!/bin/bash

for mass in 60 100 125 140 160 400 600 1000 1400 2000 2300 2900 3500; 
do

    for hdampmode in "" 1 2; 
    do 
    
        if (( $(echo $mass'<'145 | bc -l) )); 
        then 
            width=0.0041
        elif (( $(echo $mass'<'605 | bc -l) )); then 
            width=0.1
        elif (( $(echo $mass'<'2005 | bc -l) )); then
            width=1.0
        else
            width=2.0
        fi
        
        ./setup.sh $mass $width $hdampmode; 
        
    done; 
done;

