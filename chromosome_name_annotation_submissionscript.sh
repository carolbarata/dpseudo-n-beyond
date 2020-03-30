#!/bin/bash
#$ -V ## pass all environment variables to the job, VERY IMPORTANT
#$ -N chr_annotation_170 ## job name
#$ -S /bin/bash ## shell where it will run this job
#$ -j y ## join error output to normal output
#$ -cwd ## Execute the job from the current working directory
#$ -q all.q ## queue name

##################################################
# Script that changes chromosome codes to 
# chromosome names
##################################################

##### Set main variables ####
ALIGNDIR=/storage/home/users/cdcbrb/dpseudo_data/new_ref_mapping/

##### Change annotation to meaningful chromosome names #####
for SAMPLE in $(find ${ALIGNDIR} -maxdepth 1 -type d -name 'G170_BL_F');
do
	echo "Reading in ${SAMPLE}"
	sed -i.bup '/CM009975.1/ s//XL/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM009975.1/ s//XL/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	sed -i.bup '/CM009976.1/ s//XR/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM009976.1/ s//XR/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	sed -i.bup '/CM009979.1/ s//X/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM009979.1/ s//X/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	sed -i.bup '/CM000070.4/ s//2/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM000070.4/ s//2/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	sed -i.bup '/CM000071.4/ s//3/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM000071.4/ s//3/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	sed -i.bup '/CM009977.1/ s//4/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM009977.1/ s//4/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	sed -i.bup '/CM009978.1/ s//5/g' ${SAMPLE}/*_bwa_both_callers_snp_data.txt
	sed -i.bup '/CM009978.1/ s//5/g' ${SAMPLE}/*_novoalign_both_callers_snp_data.txt
	echo "${SAMPLE} done"
#	sed -i 's/ADE02000016.1/contig16/g' ${SAMPLE}/*_qual_data.txt
#	sed -i 's/ADE02000017.1/contig17/g' ${SAMPLE}/*_qual_data.txt
#	sed -i 's/ADE02000018.1/contig18/g' ${SAMPLE}/*_qual_data.txt

done
