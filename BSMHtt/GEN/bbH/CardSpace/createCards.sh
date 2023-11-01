#!/bin/bash

for mass in 100 125 1000;
# for mass in 50 60 70 80 90 95 100 105 110 120 125 130 140 160 180 200 250 300 400 500 600 800 900 1000 1100 1200 1400 1600 1800 2000 2500 3000 3500;
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

