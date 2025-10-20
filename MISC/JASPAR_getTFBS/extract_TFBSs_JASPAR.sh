#!/bin/bash

usage() {
  echo ""
  echo "Usage: $(basename $0) -i INPUT BED -b BIGBED FILE [-o OUTPUT] [-t TFs] [-m MATRIX IDs] [-s TFBS SCORE THRESHOLD]"
  echo ""
  echo "Intersect a BED file with a JASPAR TFBS bigBed file. Optionally, provide a list of TF gene symbols, JASPAR matrix IDs or TFBS score threshold to filter the results."
  echo ""
  echo "Arguments:"
  echo ""
  echo "    -i INPUT BED: path to the input BED file."
  echo "    -b BIGBED FILE: path to the bigBed file."
  echo "    -o OUTPUT: path to the output file. Results will be sent to stdout if no output is given."
  echo "    -t TFs: file containing a list of TF gene symbols of interest (case-sensitive) separated by a new line."
  echo "    -m MATRIX IDs: file containing a list of JASPAR matrix IDs of interest separated by a new line."
  echo "    -s TFBS SCORE THRESHOLD: TFBS score threshold."
  echo "    -p PROCESSORS: number of cores to run in parallel (default = 2)."
  echo "    -h HELP: show this help message."
  echo ""
  exit 0
}

#bedtools=./bedtools


#Get full path of the script - useful to make sure the script can be launched from anywhere
DIR=$( realpath "$( dirname "$0" )" )

TFBS_score_thr=0
output=""
num_cores=2

PASSED_ARGS="$@"

if [ ! -z "${PASSED_ARGS}" ]
then

  while getopts "i:b:o:t:m:s:p:h" options; do
    case "${options}" in
      i)
        input_regions=${OPTARG}
        if [ ! -e $input_regions ]
        then
          echo "Input BED file does not exist! Exiting."
          exit 1
        fi
        ;;
      b)
        bigBed_file=${OPTARG}
        if [ ! -e $bigBed_file ]
        then
          echo "Input bigBed file does not exist! Exiting."
          exit 1
        fi
        ;;
      o)
	      output=${OPTARG}
      	if [ -e $output ]
      	then
      	  echo "Output file already exists! Exiting."
	        exit 1
      	fi
	      ;;
      t)
        TFs=${OPTARG}
        if [ ! -e $TFs ]
        then
          echo "Input TFs file does not exist! Exiting."
          exit 1
        else
          TFs_comma_separated=$(cat $TFs | paste -s -d, -)
          TFs_array=(${TFs_comma_separated//","/" "}) #Separate the comma-separated input vector and put it into an array
        fi
        ;;
      m)
        matrix_IDs=${OPTARG}
        if [ ! -e $matrix_IDs ]
        then
          echo "Input matrix ID file does not exist! Exiting."
          exit 1
        else
          matrix_IDs_comma_separated=$(cat $matrix_IDs | paste -s -d, -)
          matrix_IDs_array=(${matrix_IDs_comma_separated//","/" "}) #Separate the comma-separated input vector and put it into an array
        fi
        ;;
      s)
        TFBS_score_thr=${OPTARG}
        ;;
      p)
        num_cores=${OPTARG}
        ;;
      h)
        usage
        ;;
      *)
        usage
        ;;
    esac
  done

else

  usage

fi

## ---- ##
## Main ##
## ---- ##

if [ ! -z $output ]
then

  output_dir=$(dirname $output)

else

  output_dir=$(mktemp -d)

fi

# Extract data and apply filtering only when needed
# A bit redundant, but efficient

# Solution to intersect BED with bigBed is from 
# https://bioinformatics.stackexchange.com/questions/18552/methods-to-work-directly-with-bigbed-file

cat $input_regions \
| parallel --will-cite --lb -j ${num_cores} --colsep '\t' \
${DIR}/process_bigBed.sh {1} {2} {3} $bigBed_file $output_dir

cat ${output_dir}/.tmp/* \
| ./bedtools sort -i stdin \
> ${output_dir}/extraction_results.bed

rm -r ${output_dir}/.tmp

if [ ${#TFs_array[@]} -gt 0 ]
then

  awk \
    -v OFS='\t' \
    -v column=7 \
    -v array=${TFs_comma_separated} \
    -f ${DIR}/filter_results.awk \
    ${output_dir}/extraction_results.bed \
  > ${output_dir}/extraction_results.bed.tmp

  mv \
    ${output_dir}/extraction_results.bed.tmp \
    ${output_dir}/extraction_results.bed

fi

if [ $TFBS_score_thr -gt 0  ]
then

  awk \
    -v OFS='\t' \
    -v column=5 \
    -v score=${TFBS_score_thr} \
    -f ${DIR}/filter_results_score.awk \
    ${output_dir}/extraction_results.bed \
  > ${output_dir}/extraction_results.bed.tmp

  mv \
    ${output_dir}/extraction_results.bed.tmp \
    ${output_dir}/extraction_results.bed

fi

if [ ${#matrix_IDs_array} -gt 0 ]
then

    awk \
    -v OFS='\t' \
    -v column=4 \
    -v array=${matrix_IDs_comma_separated} \
    -f ${DIR}/filter_results.awk \
    ${output_dir}/extraction_results.bed \
  > ${output_dir}/extraction_results.bed.tmp

  mv \
    ${output_dir}/extraction_results.bed.tmp \
    ${output_dir}/extraction_results.bed

fi

if [ ! -z $output ]
then

  mv \
    ${output_dir}/extraction_results.bed \
    $output

else

  cat ${output_dir}/extraction_results.bed \
  >&1

fi


