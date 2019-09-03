#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Pipeline that gets the number of initial reads
# and also those that passed Trimmomatic quality
# filters and produces a single output file for
# all available forward paired trimmed reads
##################################################

USAGE="$(basename "$0") [-h] [-r -t -o] --

where:\n
    -h show this help text\n
    -r full path to forward raw read fastq file\n
    -t full path to trimmed fastq file\n
    -s sample name\n
    -o full path to output file (e.g. filename.counttrimmed)"


while (( "$#" )); do
case "$1" in
    -r|--rawreadfile)
    RAWREADFILE="$2"
    shift 2
    ;;
    -t|--trimmedreadfile)
    TRIMMEDFILE="$2"
    shift 2
    ;;
    -s|--samplename)
    SAMPLE_NAME="$2"
    shift 2
    ;;
    -o|--outfile)
    OUTFILE="$2"
    shift 2
    ;;
    -h|--help)
    echo -e $USAGE
    shift 2
    ;;
    --) # end argument parsing
	shift
	break
	;;
    *)
    ;;
esac
done

##### Get final number of reads from forward trimmed read file
TRIMMED_READS=$(zcat ${TRIMMEDFILE} | awk '{s++}END{print s/4}')

##### Get initial number of reads from fq file
RAW_READS=$(zcat ${RAWREADFILE} | awk '{s++}END{print s/4}')

##### Add value + sample name to file
echo "${SAMPLE_NAME},${RAW_READS},${TRIMMED_READS}" > ${OUTFILE}

