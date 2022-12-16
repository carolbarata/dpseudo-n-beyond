#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N filt_matchmap_202 ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q centos7.q ## queue name

##################################################
# Script that filters SNP data files and matches
# variants called by both bwa mem and novoalign
##################################################

#    #						   #
#     #						   #
# ###### REMEMBER TO CHANGE TIME POINT PARAM '-t 1'#
#     #						   #
#    #						   #

conda activate bcftools_env

SCRIPT_DIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/vcf_analysis/
MAINDIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/G202*_F/
PYPATH=/shelf/apps/cdcbrb/conda/envs/bcftools_env/bin/python3.5

for FILE in $(find ${MAINDIR} -type f -name 'G*_F_bwa_both_callers_snp_data_XL.txt');
do
	FILE_PREFIX=${FILE%\_bwa_both_callers_snp_data_XL.txt}
	echo ${FILE_PREFIX}
	${PYPATH} ${SCRIPT_DIR}main_sample_snp_filtering_pipeline.py -b ${FILE} -n ${FILE_PREFIX}_novoalign_both_callers_snp_data_XL.txt -t 5 -s 40 -o ${FILE_PREFIX}_filtered_snp_data_XL.txt
	${PYPATH} ${SCRIPT_DIR}main_sample_snp_filtering_pipeline.py -b ${FILE_PREFIX}_bwa_both_callers_snp_data_XR.txt -n ${FILE_PREFIX}_novoalign_both_callers_snp_data_XR.txt -t 5 -s 40 -o ${FILE_PREFIX}_filtered_snp_data_XR.txt
	${PYPATH} ${SCRIPT_DIR}main_sample_snp_filtering_pipeline.py -b ${FILE_PREFIX}_bwa_both_callers_snp_data_2.txt -n ${FILE_PREFIX}_novoalign_both_callers_snp_data_2.txt -t 5 -s 40 -o ${FILE_PREFIX}_filtered_snp_data_2.txt
	${PYPATH} ${SCRIPT_DIR}main_sample_snp_filtering_pipeline.py -b ${FILE_PREFIX}_bwa_both_callers_snp_data_3.txt -n ${FILE_PREFIX}_novoalign_both_callers_snp_data_3.txt -t 5 -s 40 -o ${FILE_PREFIX}_filtered_snp_data_3.txt
	${PYPATH} ${SCRIPT_DIR}main_sample_snp_filtering_pipeline.py -b ${FILE_PREFIX}_bwa_both_callers_snp_data_4.txt -n ${FILE_PREFIX}_novoalign_both_callers_snp_data_4.txt -t 5 -s 40 -o ${FILE_PREFIX}_filtered_snp_data_4.txt
	${PYPATH} ${SCRIPT_DIR}main_sample_snp_filtering_pipeline.py -b ${FILE_PREFIX}_bwa_both_callers_snp_data_5.txt -n ${FILE_PREFIX}_novoalign_both_callers_snp_data_5.txt -t 5 -s 40 -o ${FILE_PREFIX}_filtered_snp_data_5.txt
done

conda deactivate
