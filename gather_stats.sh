#!/bin/bash
SAMPLE=$1
find $1 -name summary.txt | xargs cat 
