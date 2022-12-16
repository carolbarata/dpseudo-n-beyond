#############################################
# Script that produces config files and runs
# Bait-ER on any sync file
#############################################

##### Set main variables ####
BAITER_DIR=/path/to/Bait-ER/
CONFIG_DIR=/path/to/Bait-ER/config/files/
SYNC_DIR=/path/to/sync/files/
OUTPUT_DIR=/path/to/output/directory/

for FILE in $(find ${SYNC_DIR} -maxdepth 1 -type f -name 'R4T3_Tfirst3_chr4_*_snp_only_nohead.sync');
do
    ##### Produce config for each subset #####
    LOCI=$(wc -l ${FILE} | awk '{print $1;}')
    CONFIG_FILE=${CONFIG_DIR}${FILE##*/}
    CONFIG_FILE=${CONFIG_FILE%.sync}.cf
    OUTPUT_FILE=${OUTPUT_DIR}${FILE##*/}
    OUTPUT_FILE=${OUTPUT_FILE%.sync}.baiter

    # Hash if Tall or Tfirst3
    #echo -e "Sync_file\t${FILE}\nHeader\t0\nColumns_order\t0\nNumber_replicates\t4\nTime_points\t21,61,114,162,200\nNumber_loci\t${LOCI}\nPopulation_size\t160\nPrior_parameters\t0.001,0.001\nOutput_file\t${OUTPUT_FILE}" >> ${CONFIG_FILE}
    echo -e "Sync_file\t${FILE}\nHeader\t0\nColumns_order\t0\nNumber_replicates\t4\nTime_points\t21,61,114\nNumber_loci\t${LOCI}\nPopulation_size\t160\nPrior_parameters\t0.001,0.001\nOutput_file\t${OUTPUT_FILE}" >> ${CONFIG_FILE}

    ##### Run Bait-ER on each subset #####
    cmd="${BAITER_DIR}baiter ${CONFIG_FILE}"
    echo ${cmd}
    eval ${cmd}
done
