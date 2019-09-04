#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Pipeline that can parse read data using
# Trimmomatic, map it using bwa mem or NovoAlign,
# and compute some read and mapping quality stats
##################################################

USAGE="$(basename "$0") [-h] [-i1 -i2 -u -dt -df -dm -dn -ds -dg -i -r -o] --

where:\n
    -h show this help text\n
    -p full path to pipeline directory\n
    -i1 full path to forward fastq read file\n
    -i2 full path to reverse fastq read file\n
    -u whether trimmed reads should be mapped to the ref (TRUE/FALSE)\n
    -dt full path to Trimmomatic directory (NA if reads are not to be trimmed)\n
    -df full path to FastQC directory (NA if reads are not to go through QC)\n
    -dm full path to bwa mem directory (NA if bwa mem is not to be used)\n
    -dn full path to NovoAlign directory (NA if NovoAlign is not to be used)\n
    -dj full path to Java 1.8 directory\n
    -ds full path to samtools directory\n
    -dg full path to GATK directory\n
    -db full path to sambamba directory\n
    -i whether the reference genome is already indexed or not (TRUE/FALSE)\n
    -r full path to reference genome\n
    -o full path to output directory (e.g. /full/path/to/outdir/)"


while (( "$#" )); do
case "$1" in
    -p|--pipelinedirectory)
    PIPEDIR="$2"
    shift 2
    ;;
    -i1|--forwardfastqfile)
    FASTQFILE1="$2"
    shift 2
    ;;
    -i2|--reversefastqfile)
    FASTQFILE2="$2"
    shift 2
    ;;
    -u|--usetrimmed)
    USETRIMMED="$2"
    shift 2
    ;;
    -dt|--trimmomaticdir)
    TRIMMOMATICDIR="$2"
    shift 2
    ;;
    -df|--fastqcdir)
    FASTQCDIR="$2"
    shift 2
    ;;
    -dm|--bwamemdir)
    BWADIR="$2"
    shift 2
    ;;
    -dn|--novoaligndir)
    NOVOALIGNDIR="$2"
    shift 2
    ;;
    -dj|--java18dir)
    JAVA18DIR="$2"
    shift 2
    ;;
    -ds|--samtoolsdir)
    SAMTOOLSDIR="$2"
    shift 2
    ;;
    -dg|--gatkdirectory)
    GATKDIR="$2"
    shift 2
    ;;
    -db|--sambambadirectory)
    SAMBAMBADIR="$2"
    shift 2
    ;;
    -i|--index)
    INDEXING="$2"
    shift 2
    ;;
    -r|--refgenome)
    REFGENOME="$2"
    shift 2
    ;;
    -o|--outdir)
    OUTDIR="$2"
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


########## Mapping pipeline ##########

## NB: should take one sample at a time ##

##### Set variables and directories
SAMPLENAME=$(echo ${FASTQFILE1} | grep -oE '[^/]+_F' | tail -n1)

if [ ${USETRIMMED} = "TRUE" ];
then 
	OUTSAMPLENAME=${SAMPLENAME}_trim
elif [ ${USETRIMMED} = "FALSE" ];
then
	OUTSAMPLENAME=${SAMPLENAME}
fi

[ -d "${OUTDIR}${SAMPLENAME}" ] || mkdir ${OUTDIR}${SAMPLENAME}

##### (0) Trim reads - OPTIONAL
if [ ${TRIMMOMATICDIR} != "NA" ];
then
	# Filters raw reads using Trimmomatic
	${PIPEDIR}trim_raw_reads_with_trimmomatic.sh \
	-d ${TRIMMOMATICDIR} \
	-i1 ${FASTQFILE1} -i2 ${FASTQFILE2} \
	-s ${OUTDIR}${SAMPLENAME}/${SAMPLENAME}
	
	#### (0.1) Compare trimmed reads with raw data
	### (0.1.1) Count number of reads before/after trimming
	# Loops over samples found in a given directory
	${PIPEDIR}get_initial_and_trimmed_read_number.sh \
	-r ${FASTQFILE1} \
	-t ${OUTDIR}${SAMPLENAME}/${SAMPLENAME}_forward_paired.fq.gz \
	-s ${SAMPLENAME} \
	-o ${OUTDIR}${SAMPLENAME}/${SAMPLENAME}.countreads
	
	### (0.1.2) Get read length distn before/after trimming
	# Loops over samples found in a given directory
	${PIPEDIR}get_initial_and_trimmed_read_lengths.sh \
	-r ${FASTQFILE1} \
	-t ${OUTDIR}${SAMPLENAME}/${SAMPLENAME}_forward_paired.fq.gz \
	-o ${OUTDIR}${SAMPLENAME}/${SAMPLENAME}.readlengths
	
	## (0.1.2.1) R-script to produce read length distn plot
	
	Rscript --vanilla ${PIPEDIR}raw_vs_trimmed_read_length_plot.R \
	${OUTDIR}${SAMPLENAME}/${SAMPLENAME}.readlengths \
	${OUTDIR}${SAMPLENAME}/${SAMPLENAME}.readlengths.png

        FASTQFILE1=${OUTDIR}${SAMPLENAME}/${SAMPLENAME}_forward_paired.fq.gz
        FASTQFILE2=${OUTDIR}${SAMPLENAME}/${SAMPLENAME}_reverse_paired.fq.gz
	[ -s "${OUTDIR}${SAMPLENAME}/${SAMPLENAME}_forward_paired.fq.gz" ] && echo "Trimming was probably successful :)" || echo "Whoops, trimming doesnt seem to have worked."
fi

##### (1) FastQC input reads - OPTIONAL
if [ ${FASTQCDIR} != "NA" ] && [ ${USETRIMMED} = "TRUE" ];
then
	${FASTQCDIR}fastqc ${OUTDIR}${SAMPLENAME}/*.fq.gz
elif [ ${FASTQCDIR} != "NA" ] && [${USETRIMMED} = "FALSE" ];
then
	${FASTQCDIR}fastqc ${OUTDIR}${OUTSAMPLENAME}/*.fq.gz
fi

##### (2) Map reads to reference genome
if [[ ${REFGENOME} == *"fna.gz" ]];
then
	UNZIPREF=$(echo ${REFGENOME} | grep -oE '[^/]+.fna')
elif [[ ${REFGENOME} == *"fasta.gz" ]];
then
	UNZIPREF=$(echo ${REFGENOME} | grep -oE '[^/]+.fasta')
fi

REFDIR=$(echo ${REFGENOME} | sed 's:/[^/]*$::')

if [ ${BWADIR} != "NA" ]
then
	if [ ${INDEXING} == "FALSE" ]
	then
		${BWADIR}bwa index ${REFGENOME}
	fi
	#${BWADIR}bwa mem ${REFGENOME} ${FASTQFILE1} ${FASTQFILE2} | gzip > "${OUTDIR}${SAMPLENAME}/${OUTSAMPLENAME}.sam.gz"
	SAMFILE=${OUTDIR}${SAMPLENAME}/${OUTSAMPLENAME}.sam.gz
	TESTSAM=$( head -c 1 < <(zcat "${SAMFILE}") )	
	[ -s "${SAMFILE}" -a "${TESTSAM}" == "@" ] && echo "Looks like some reads got mapped okay." || echo "Mapping seems to have gone wrong, oops."
elif [ ${NOVOALIGNDIR} != "NA" ]
then
	if [ ${INDEXING} == "FALSE" ]
	then
		${NOVOALIGNDIR}novoindex ${REFDIR}/novoalignref.nix ${REFDIR}/${UNZIPREF}
	fi
	${NOVOALIGNDIR}novoalign \
	-d ${REFDIR}/novoalignref.nix \
	-f ${FASTQFILE1} ${FASTQFILE2} \
	-o SAM > "${OUTDIR}${SAMPLENAME}/${OUTSAMPLENAME}_novoalign.sam"
	gzip "${OUTDIR}${SAMPLENAME}/${OUTSAMPLENAME}_novoalign.sam"
	SAMFILE="${OUTDIR}${SAMPLENAME}/${OUTSAMPLENAME}_novoalign.sam.gz"
	TESTSAM=$( head -c 1 < <(zcat "${SAMFILE}") )
	[ -s "${SAMFILE}" -a "${TESTSAM}" == "@" ] && echo "Looks like some reads got mapped okay." || echo "Mapping seems to have gone wrong, oops."
fi

#### (3) Realign around indels and compute mapping stats
# Parses a sam file
# realigns around indels using GATK
# and computes mapping stats using samtools
SAMSAMPLENAME=$(echo ${SAMFILE%\.sam.gz} | grep -oE '[^/]+' | tail -n1)
${PIPEDIR}parsing_mapped_reads_with_samtools.sh \
-d ${SAMTOOLSDIR} \
-g ${GATKDIR} \
-j ${JAVA18DIR} \
-b ${SAMBAMBADIR} \
-s ${SAMSAMPLENAME} \
-i ${SAMFILE} \
-r ${REFDIR}/${UNZIPREF} \
-u ${USETRIMMED} \
-o ${OUTDIR}

# Extracts mapping quality score using samtools
${PIPEDIR}mapq_from_samfiles.sh \
-d ${SAMTOOLSDIR} \
-i ${SAMFILE} \
-o ${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.mapq

# R-script to plot coverage (aka sequencing depth)
Rscript --vanilla ${PIPEDIR}coverage_distn_plot.R \
${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.nodup.sorted.bam.depth \
${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.depth.png

        [ -s "${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.mapq" -a -s "${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.depth.png" ] && echo "Wonderful, parsing seems to have gone just right." || echo "Hmm, something seems to have gone wrong. Please check why your mapped reads havent been parsed correctly."


if [ ${USETRIMMED} = "TRUE" ] && [ ${BWADIR} != "NA" ] || [ ${NOVOALIGNDIR} != "NA" ];
then
	#### (3.1) Compare with raw data mapping stats - OPTIONAL
	### (3.1.1) Compare number of mapped reads before/after trimming
	RAW_MAPPED=$(${SAMTOOLSDIR}samtools view -c -F 4 ${SAMFILE})
	TRIM_MAPPED=$(${SAMTOOLSDIR}samtools view -c -F 4 ${SAMFILE})
	echo "${SAMPLENAME},${RAW_MAPPED},${TRIM_MAPPED}" > ${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.mappedreads
	### (3.1.2) Compare mapping quality scores before/after trimming
        Rscript --vanilla ${PIPEDIR}raw_vs_trimmed_mapped_read_distn_plot.R \
        ${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.mapq \
        ${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.mapq \
        ${OUTDIR}${SAMPLENAME}/${SAMSAMPLENAME}.mapq.png
fi
