#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Script that uses bcftools or freebayes to
# perform variant calling
# NB(1): input bamfile should refer to a single
# treatment (all replicates, all time points)
# NB(2): if freebayes is caller of choice,
# script finds common variants between
# bcftools&fb datasets and keeps bcftools entries
##################################################

USAGE="$(basename "$0") [-h] [-id -b -s -fb -f -o -r -i] --

where:\n
    -h show this help text\n
    -id sample ID/name\n
    -m name of program used for read mapper (bwa or novoalign)\n
    -b full path to bcftools directory\n
    -s full path to samtools directory\n
    -z full path to bgzip directory\n
    -fb full path to freebayes directory\
    -f full path to input bam file\n
    -o full path to output directory\n
    -r full path to reference fasta file\n
    -i whether the reference genome is already indexed or not (TRUE or FALSE)"


while (( "$#" )); do
case "$1" in
    -id|--sampleid)
    SAMPLEID="$2"
    shift 2
    ;;
    -m|--mapper)
    MAPPER="$2"
    shift 2
    ;;
    -b|--bcftoolsdirectory)
    BCFTOOLSDIR="$2"
    shift 2
    ;;
    -s|--samtoolsdirectory)
    SAMTOOLSDIR="$2"
    shift 2
    ;;
    -z|--bgzipdirectory)
    BGZIPDIR="$2"
    shift 2
    ;;
    -fb|--freebayesdirectory)
    FREEBAYESDIR="$2"
    shift 2
    ;;
    -f|--bamfile)
    BAMFILE="$2"
    shift 2
    ;;
    -o|--outputdirectory)
    OUTDIR="$2"
    shift 2
    ;;
    -r|--refgenome)
    REFGEN="$2"
    shift 2
    ;;
    -i|--indexing)
    INDEXING="$2"
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


########## SNP calling pipeline ##########

## NB: should take one time point/replicate at a time ##

## Index reference genome
if [ ${INDEXING} = "FALSE" ]; then ${SAMTOOLSDIR}samtools faidx ${REFGEN}; fi

## Rename SAMPLEID
SAMPLEID=${SAMPLEID}_${MAPPER}

##### SNP calling #####
if [ ${BCFTOOLSDIR} != "NA" ] && [ ${FREEBAYESDIR} = "NA" ];
then
	${BCFTOOLSDIR}bcftools mpileup -Ou \
        -f ${REFGEN} ${BAMFILE} \
	-a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR |
       ${BCFTOOLSDIR}bcftools call -m -Ob \
       -o ${OUTDIR}${SAMPLEID}_unfiltered.vcf.gz
	[ -s "${OUTDIR}${SAMPLEID}_unfiltered.vcf.gz" ]  || echo "Hmm, there's no output VCF file, might want to check what went wrong with SNP calling."
	TESTVCF=$(${BCFTOOLSDIR}bcftools view -h ${OUTDIR}${SAMPLEID}_unfiltered.vcf.gz | head -n1 | grep -oE '^.*VCF')
	[ ${TESTVCF} = "##fileformat=VCF" ] && echo "Woohoo, SNP calling seems to have been successful."
	# Index .vcf.gz file
	${BCFTOOLSDIR}bcftools index -f ${OUTDIR}${SAMPLEID}_unfiltered.vcf.gz
	# Compute vcf stats with bcftools
	${BCFTOOLSDIR}bcftools stats ${OUTDIR}${SAMPLEID}_unfiltered.vcf.gz > ${OUTDIR}/${SAMPLEID}_unfiltered.vcfstats

	## Filter SNPs based on quality criteria
	# Min base quality 30
	${BCFTOOLSDIR}bcftools view \
	        -i '%QUAL>=30' \
	        -U \
	        -Oz ${OUTDIR}${SAMPLEID}_unfiltered.vcf.gz \
	         > ${OUTDIR}${SAMPLEID}_filtered.vcf.gz
	
	[ -s "${OUTDIR}${SAMPLEID}_filtered.vcf.gz" ] || echo "Whoops, there's no final filtered VCF file, SNP filtering didnt work."
	# Compute vcf stats with bcftools
        ${BCFTOOLSDIR}bcftools stats ${OUTDIR}${SAMPLEID}_filtered.vcf.gz > ${OUTDIR}/${SAMPLEID}_filtered.vcfstats
	TESTFILT=$(${BCFTOOLSDIR}bcftools view -h ${OUTDIR}${SAMPLEID}_filtered.vcf.gz | head -n1 | grep -oE '^.*VCF')
	[ ${TESTFILT} = "##fileformat=VCF" ] && echo "Phew, bcftools SNP filtering is finally over. Congratulations."
fi

if [ ${FREEBAYESDIR} != "NA" ];
then
	${FREEBAYESDIR}freebayes-v1.3.1 \
	-f ${REFGEN} \
	--report-monomorphic \
	--pooled-continuous  \
	${BAMFILE} > ${OUTDIR}${SAMPLEID}_freebayes.vcf
	${BGZIPDIR}bgzip ${OUTDIR}${SAMPLEID}_freebayes.vcf
	[ -s "${OUTDIR}${SAMPLEID}_freebayes.vcf.gz" ] || echo "Hmm, there's no output VCF file, might want to check what went wrong with SNP calling."
        ${BCFTOOLSDIR}bcftools index -f ${OUTDIR}${SAMPLEID}_freebayes.vcf.gz
	# Compute vcf stats with bcftools
        ${BCFTOOLSDIR}bcftools stats ${OUTDIR}${SAMPLEID}_freebayes.vcf.gz > ${OUTDIR}${SAMPLEID}_freebayes.vcfstats
	TESTFBVCF=$(${BCFTOOLSDIR}bcftools view -h ${OUTDIR}${SAMPLEID}_freebayes.vcf.gz | head -n1 | grep -oE '^.*VCF')
	[ ${TESTFBVCF} = "##fileformat=VCF" ] && echo "Woohoo, freebayes SNP calling seems to have been successful."
	${BCFTOOLSDIR}bcftools index -f ${OUTDIR}${SAMPLEID}_filtered.vcf.gz
	${BCFTOOLSDIR}bcftools isec -p ${OUTDIR%\/} ${OUTDIR}${SAMPLEID}_filtered.vcf.gz ${OUTDIR}${SAMPLEID}_freebayes.vcf.gz
	[ -s "${OUTDIR}0002.vcf" ] && echo "VCF intersection seems to have gone through but please make sure your files arent empty." || echo "Whoops, intersection didnt work."
	mv ${OUTDIR}0002.vcf ${OUTDIR}${SAMPLEID}_both_callers.vcf
	mv ${OUTDIR}README.txt ${OUTDIR}README_${MAPPER}_both_callers.txt
	rm ${OUTDIR}0000.vcf ${OUTDIR}0001.vcf ${OUTDIR}0003.vcf
	${BGZIPDIR}bgzip ${OUTDIR}${SAMPLEID}_both_callers.vcf
	${BCFTOOLSDIR}bcftools index -f ${OUTDIR}${SAMPLEID}_both_callers.vcf.gz
	${BCFTOOLSDIR}bcftools stats ${OUTDIR}${SAMPLEID}_both_callers.vcf.gz > ${OUTDIR}${SAMPLEID}_both_callers.vcfstats
fi

