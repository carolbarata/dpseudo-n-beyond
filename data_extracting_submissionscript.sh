#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N snp_extract_g ## job name
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
ALIGNDIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/
SCRIPTDIR=/storage/home/users/cdcbrb/dpseudo_mapping_pipeline/

##### Count variants with bcftools #####
touch ${ALIGNDIR}snp_count_per_sample_bwa_both_callers.counts
touch ${ALIGNDIR}snp_count_per_sample_novoalign_both_callers.counts
for SAMPLE in $(find ${ALIGNDIR} -maxdepth 1 -type d -name 'G*_F');
do
	SAMPLENAME=$(echo ${SAMPLE} | grep -oE '[^/]+_F')
	BWAVCF=${SAMPLE}/*_bwa_both_callers.vcf.gz
	NOVOALIGNVCF=${SAMPLE}/*_novoalign_both_callers.vcf.gz
	${SCRIPTDIR}main_data_extracting_pipeline.sh -id ${SAMPLENAME}_bwa_both_callers -i ${BWAVCF} -b ${BCFTOOLSDIR} -z ${BGZIPDIR} -o ${SAMPLE} -c ${ALIGNDIR}snp_count_per_sample_bwa_both_callers.counts
	${SCRIPTDIR}main_data_extracting_pipeline.sh -id ${SAMPLENAME}_novoalign_both_callers -i ${NOVOALIGNVCF} -b ${BCFTOOLSDIR} -z ${BGZIPDIR} -o ${SAMPLE} -c ${ALIGNDIR}snp_count_per_sample_novoalign_both_callers.counts

done

conda deactivate
