#!/bin/bash
set -e
set -u
set -o pipefail

##################################################
# Pipeline for computing mapping quality
# control and other mapping stats with samtools
##################################################

USAGE="$(basename "$0") [-h] [ -d -g -j -b -s -i -r -o] --

where:\n
    -h show this help text\n
    -d full path to samtools directory\n
    -g full path to gatk directory\n
    -j full path to Java 1.8 directory\n
    -b full path to sambamba directory\n
    -s sample name\n
    -i full path to samfile input file\n
    -r full path to reference genome file\n
    -o full path to output directory"


while (( "$#" )); do
case "$1" in
    -d|--samtoolsdirectory)
    SAMTOOLSDIR="$2"
    shift 2
    ;;
    -g|--gatkdirectory)
    GATKDIR="$2"
    shift 2
    ;;
    -j|--javadirectory)
    JAVA18DIR="$2"
    shift 2
    ;;
    -b|--sambambadirectory)
    SAMBAMBADIR="$2"
    shift 2
    ;;
    -s|--samplename)
    SAMPLENAME="$2"
    shift 2
    ;;
    -i|--infile)
    SAMFILE="$2"
    shift 2
    ;;
    -r|--refgenome)
    REFGENOME="$2"
    shift 2
    ;;
    -o|--outdirectory)
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

##### Set main variables #####
FINBAMEXT=".nodup.sorted.realign.bam"
IDXEXT=".nodup.sorted.bam.idxstats"
STATSEXT=".nodup.sorted.bam.stats"
FLAGEXT=".nodup.sorted.bam.flagstat"
COVEXT=".nodup.sorted.bam.depth"

BAMFILE=${OUTDIR}/${SAMPLENAME}_temp.bam
MATE_BAMFILE=${OUTDIR}/${SAMPLENAME}_temp.fixmate.bam
FIX_BAMFILE=${OUTDIR}/${SAMPLENAME}_temp.nounmap.bam
NODUP_BAMFILE=${OUTDIR}/${SAMPLENAME}_temp.nodup.bam
SORTED_BAMFILE=${OUTDIR}/${SAMPLENAME}_temp.sorted.bam
GROUP_BAMFILE=${OUTDIR}/${SAMPLENAME}_temp.groupinfo.bam
INTERVALS_FILE=${OUTDIR}/${SAMPLENAME}_temp.intervals
REALBAMFILE=${OUTDIR}/${SAMPLENAME}_temp.realign.bam
FINBAMFILE=${OUTDIR}/${SAMPLENAME}${FINBAMEXT}
IDXSTATS_FILE=${OUTDIR}/${SAMPLENAME}${IDXEXT}
STATS_FILE=${OUTDIR}/${SAMPLENAME}${STATSEXT}
FLAGSTAT_FILE=${OUTDIR}/${SAMPLENAME}${FLAGEXT}
COV_FILE=${OUTDIR}/${SAMPLENAME}${COVEXT}

##### Compute alignment stats #####
cd ${SAMTOOLSDIR}

##### Set up java 1.8 directory #####
export JAVA_HOME=${JAVA18DIR}
export PATH="$JAVA_HOME:$PATH"

##### Produce bam file from sam alignment #####
./samtools view -S -b -q 40 $SAMFILE > $BAMFILE

##### Fix mate pairs with picard #####
${JAVA18DIR}java -jar ${GATKDIR}picard.jar FixMateInformation \
I=${BAMFILE} \
O=${MATE_BAMFILE}

##### Fix mate pairs with samtools #####
./samtools fixmate -rcm $MATE_BAMFILE $FIX_BAMFILE # fill in mate coordinates, ISIZE and mate related flags

##### Remove duplicate reads with sambamba #####
${SAMBAMBADIR}sambamba markdup -r ${MATE_BAMFILE} ${NODUP_BAMFILE}

./samtools sort $NODUP_BAMFILE > $SORTED_BAMFILE
./samtools index $SORTED_BAMFILE

##### Delete intermediate bam files #####
rm $BAMFILE $MATE_BAMFILE $FIX_BAMFILE $NODUP_BAMFILE

##### Add random read groups to ${BAMFILE} #####
${JAVA18DIR}java -jar ${GATKDIR}picard.jar AddOrReplaceReadGroups \
I=${SORTED_BAMFILE} \
O=${GROUP_BAMFILE} \
RGID=${SAMPLENAME} \
RGLB=L1 \
RGPL=illumina \
RGPU=NONE \
RGSM=S1

./samtools index $GROUP_BAMFILE

##### Create target intervals file #####
${JAVA18DIR}java -jar ${GATKDIR}GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R ${REFGENOME} \
-I ${GROUP_BAMFILE} \
-o ${INTERVALS_FILE}

##### Realign around indels #####
${JAVA18DIR}java -jar ${GATKDIR}GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ${REFGENOME} \
-I ${GROUP_BAMFILE} \
-targetIntervals ${INTERVALS_FILE} \
-o ${REALBAMFILE}

./samtools sort $REALBAMFILE > $FINBAMFILE
./samtools index $FINBAMFILE # index a coordinate-sorted BAM file for fast random access
./samtools idxstats $FINBAMFILE > $IDXSTATS_FILE # retrieve and print stats in the index file
./samtools stats $FINBAMFILE > $STATS_FILE
./samtools flagstat $FINBAMFILE > $FLAGSTAT_FILE
./samtools depth $FINBAMFILE > $COV_FILE # compute coverage

##### Delete intermediate files #####
rm $SORTED_BAMFILE ${SORTED_BAMFILE}.bai $GROUP_BAMFILE ${GROUP_BAMFILE}.bai $INTERVALS_FILE $REALBAMFILE
