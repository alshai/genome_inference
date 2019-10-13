#!/bin/bash
SAMPLE=$1
cat data/test_samples.txt | find $1 -name summary.txt | xargs cat 
