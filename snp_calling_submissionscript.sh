#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N novo_fb_21m ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q centos7.q ## queue name

##################################################
# Script that calls snps using bcftools
##################################################

conda activate bcftools_env

##### Set main variables ####
NEWREF=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/dpse-xchr-ucberk-1.0.fasta
#OLDREF=
PIPELINEDIR=/storage/home/users/cdcbrb/dpseudo_mapping_pipeline/
BCFTOOLSDIR=/shelf/apps/cdcbrb/conda/envs/bcftools_env/bin/
SAMTOOLSDIR=/shelf/apps/cdcbrb/conda/envs/bcftools_env/bin/
BGZIPDIR=/shelf/apps/cdcbrb/conda/envs/bcftools_env/bin/
FREEBAYESDIR=/storage/home/users/cdcbrb/software/
ALIGNDIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/xchr/G21_M*

##### Calls variants with bcftools #####

for SAMPLE in $(find ${ALIGNDIR} -maxdepth 1 -type d);
do
	SAMPLENAME=$(echo ${SAMPLE} | grep -oE '[^/]+_F' | tail -n1)
	${PIPELINEDIR}main_snp_calling_pipeline.sh \
	-id ${SAMPLENAME} \
	-m novoalign \
	-b NA \
	-s ${SAMTOOLSDIR} \
	-z ${BGZIPDIR} \
	-fb ${FREEBAYESDIR} \
	-f ${SAMPLE}/*_trim_novoalign.nodup.sorted.realign.bam \
	-o ${SAMPLE}/ \
	-r ${NEWREF} \
	-i TRUE
done

conda deactivate
