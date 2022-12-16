#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N snp_subsetting_200p3 # job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q centos7.q ## queue name

##################################################
# Script that filters VCF files
##################################################

MAINDIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/

for SAMPLE in $(find ${MAINDIR} -maxdepth 1 -type d -name 'G200_P3_F');
do
	for FILE in $(find ${SAMPLE} -type f -name '*_both_callers_snp_data.txt');
	do
		head -n1 ${FILE} >> ${FILE%\.txt}_XL.txt && grep -E '^XL' ${FILE} >> ${FILE%\.txt}_XL.txt
		head -n1 ${FILE} >> ${FILE%\.txt}_XR.txt && grep -E '^XR' ${FILE} >> ${FILE%\.txt}_XR.txt
		head -n1 ${FILE} >> ${FILE%\.txt}_2.txt && grep -E '^2' ${FILE} >> ${FILE%\.txt}_2.txt
		head -n1 ${FILE} >> ${FILE%\.txt}_3.txt && grep -E '^3' ${FILE} >> ${FILE%\.txt}_3.txt
		head -n1 ${FILE} >> ${FILE%\.txt}_4.txt && grep -E '^4' ${FILE} >> ${FILE%\.txt}_4.txt
		head -n1 ${FILE} >> ${FILE%\.txt}_5.txt && grep -E '^5' ${FILE} >> ${FILE%\.txt}_5.txt
	done
	echo "${SAMPLE} is done"
done


