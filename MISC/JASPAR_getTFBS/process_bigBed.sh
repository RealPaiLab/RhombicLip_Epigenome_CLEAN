#!/bin/bash

chr=$1
start=$2
end=$3
bigBed_file=$4
output_dir=$5

BIGBED2BED=./kent_utilities/bigBedToBed

#Get full path of the script - useful to make sure the 
# script can be launched from anywhere
DIR=$( realpath "$( dirname "$0" )" )

mkdir -p ${output_dir}/.tmp
tmp=$(mktemp -p ${output_dir}/.tmp)

$BIGBED2BED -chrom=${chr} -start=${start} -end=${end} $bigBed_file $tmp
