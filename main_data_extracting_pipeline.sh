#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Script that uses bcftools to output a
# multi-sample pileup file
# NB: input bamfile should refer to a single
# treatment (all replicates, all time points)
##################################################

USAGE="$(basename "$0") [-h] [-id -b -s -fb -f -o -r -i] --

where:\n
    -h show this help text\n
    -id sample ID/name\n
    -i full path to all input VCF files\n
    -b full path to bcftools directory\n
    -z full path to bgzip directory\n
    -o full path to output directory\n
    -c full path to file containing SNP counts per sample"


while (( "$#" )); do
case "$1" in
    -id|--sampleid)
    SAMPLEID="$2"
    shift 2
    ;;
    -i|--inputvcf)
    VCF="$2"
    shift 2
    ;;
    -b|--bcftoolsdirectory)
    BCFTOOLSDIR="$2"
    shift 2
    ;;
    -z|--bgzipdirectory)
    BGZIPDIR="$2"
    shift 2
    ;;
    -o|--outputdirectory)
    OUTDIR="$2"
    shift 2
    ;;
    -c|--countsfile)
    COUNTSFILE="$2"
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


########## SNP VCF data extracting pipeline ##########

## NB: should take one sample/time point at a time  ##

##### Extract ALT allele variant calling & genotype calling qual scores #####
echo -e 'CHROM\tPOS\tREF\tALT\tQUAL\tGT\tDP\tADF\tADR\tAD' > ${OUTDIR}/${SAMPLEID}_snp_data.txt
#${BCFTOOLSDIR}bcftools filter --exclude TYPE="indel" ${VCF} > ${VCF}.temp
${BCFTOOLSDIR}bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%ADF\t%ADR\t%AD]\n' ${VCF} >> ${OUTDIR}/${SAMPLEID}_snp_data.txt
rm ${VCF}.temp

[ -s "${OUTDIR}/${SAMPLEID}.vcfstats" ] && [ -s "${OUTDIR}/${SAMPLEID}_snp_data.txt" ] && echo "Stats were computed and data also seems to have been extracted, woohoo!" || echo "Whoops, something's gone wrong with data extracting and/or stats computing."

##### Extract SNP number from stats file #####

echo -e "$(echo ${SAMPLEID} | grep -oE '[^/]+_F' )\t$(grep -E "^SN.*SNPs:" ${OUTDIR}/${SAMPLEID}.vcfstats | awk -F: '{print $NF}' | sed -e 's/^[ \t]*//')" >> ${COUNTSFILE}
