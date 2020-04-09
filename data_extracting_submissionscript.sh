#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N snp_extract_112 ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q centos7.q ## queue name

##################################################
# Script that extracts snp data using bcftools
##################################################

conda activate bcftools_env

##### Set main variables ####
BCFTOOLSDIR=/shelf/apps/cdcbrb/conda/envs/bcftools_env/bin/
BGZIPDIR=/shelf/apps/cdcbrb/conda/envs/bcftools_env/bin/
ALIGNDIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/xchr/
SCRIPTDIR=/storage/home/users/cdcbrb/dpseudo_mapping_pipeline/

##### Count variants with bcftools #####

ASSEMBLYLV="xchr"

if [ ${ASSEMBLYLV} = "xchr" ];
then
	BWACOUNTSFILE=${ALIGNDIR}snp_count_per_sample_bwa_both_callers_xchr.counts
	NOVOCOUNTSFILE=${ALIGNDIR}snp_count_per_sample_novoalignn_both_callers_xchr.counts
else
	BWACOUNTSFILE=${ALIGNDIR}snp_count_per_sample_bwa_both_callers.counts
	NOVOCOUNTSFILE=${ALIGNDIR}snp_count_per_sample_novoalign_both_callers.counts
fi

#### Hash two lines below if the 
#### sample you're about to parse ####
#### isn't the first ever ####

#touch ${BWACOUNTSFILE}
#touch ${NOVOCOUNTSFILE}

for SAMPLE in $(find ${ALIGNDIR} -maxdepth 1 -type d -name 'G112*_F');
do
	SAMPLENAME=$(echo ${SAMPLE} | grep -oE '[^/]+_F')
	BWAVCF=${SAMPLE}/*_bwa_both_callers.vcf.gz
	NOVOALIGNVCF=${SAMPLE}/*_novoalign_both_callers.vcf.gz
	${SCRIPTDIR}main_data_extracting_pipeline.sh -id ${SAMPLENAME}_bwa_both_callers -i ${BWAVCF} -b ${BCFTOOLSDIR} -z ${BGZIPDIR} -o ${SAMPLE} -c ${BWACOUNTSFILE} -l ${ASSEMBLYLV}
	${SCRIPTDIR}main_data_extracting_pipeline.sh -id ${SAMPLENAME}_novoalign_both_callers -i ${NOVOALIGNVCF} -b ${BCFTOOLSDIR} -z ${BGZIPDIR} -o ${SAMPLE} -c ${NOVOCOUNTSFILE} -l ${ASSEMBLYLV}

done

conda deactivate
