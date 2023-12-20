#!/bin/bash

# $1 = sorted bam file
# $2 = reference gtf file
# $3 = output file name

stringtie $1 -e -G $2 -o $3 -p 15
