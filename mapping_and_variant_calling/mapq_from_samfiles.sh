#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Pipeline that gets mapping quality score (MAPQ)
# from .sam files
##################################################

USAGE="$(basename "$0") [-h] [-d -i -o] --

where:\n
    -h show this help text\n
    -d full path to samtools directory\n
    -i full path to samfile input file\n
    -o full path to output file"


while (( "$#" )); do
case "$1" in
    -d|--samtoolsdirectory)
    SAMTOOLSDIR="$2"
    shift 2
    ;;
    -i|--infile)
    SAMFILE="$2"
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

cd ${SAMTOOLSDIR}
./samtools view ${SAMFILE} | awk -F "\t" '{print $5}' > ${OUTFILE}

