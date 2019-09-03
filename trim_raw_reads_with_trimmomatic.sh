#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Pipeline for trimming raw reads with
# Trimmomatic (for paired end reads)
##################################################

USAGE="$(basename "$0") [-h] [-d -i1 -i2 -s] --

where:\n
    -h show this help text\n
    -d full path to trimmomatic directory\n
    -i1 full path to fastqfile 1\n
    -i2 full path to fastqfile 2\n
    -s sample name"


while (( "$#" )); do
case "$1" in
    -d|--trimmomaticdirectory)
    TRIMDIR="$2"
    shift 2
    ;;
    -i1|--fastqfile1)
    FASTQFILE1="$2"
    shift 2
    ;;
    -i2|--fastqfile2)
    FASTQFILE2="$2"
    shift 2
    ;;
    -s|--samplename)
    SAMPLENAME="$2"
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

java -jar ${TRIMDIR}trimmomatic-0.38.jar PE -phred33 \
	${FASTQFILE1} \
	${FASTQFILE2} \
	${SAMPLENAME}_forward_paired.fq.gz \
	${SAMPLENAME}_forward_unpaired.fq.gz \
	${SAMPLENAME}_reverse_paired.fq.gz \
	${SAMPLENAME}_reverse_unpaired.fq.gz \
	ILLUMINACLIP:${TRIMDIR}adapters/TruSeq3-PE.fa:2:30:10 \
	LEADING:3 \
	TRAILING:3 \
	SLIDINGWINDOW:4:15 \
	MINLEN:36

