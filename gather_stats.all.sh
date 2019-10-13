#!/bin/bash

while read sample;
do 
    ./gather_stats.sh $sample
done < data/test_samples.txt
