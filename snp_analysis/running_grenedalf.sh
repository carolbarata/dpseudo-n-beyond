#############################################
# Script that runs grenedalf to produce pi  #
# estimates and Tajima's D on a sync file   #
#############################################

##### Set main variables ####
SOFTWAREDIR=/path/to/grenedalf/bin/
DATADIR=/path/to/sync/data/directory/
SAMPLE_NAMES=${DATADIR}sample_name_list_R4T5.txt
OUTPUT_DIR=/path/to/output/directory/

##### Compute pi #####
for FILE in $(find ${DATADIR} -maxdepth 1 -type f -name 'R4T5_*_snp_only_nohead.sync'); do
  FILENAME=${FILE#${DATADIR}}
  ${SOFTWAREDIR}grenedalf diversity --sync-file ${DATADIR}${FILENAME} \
                        --window-width 250000 \
                        --window-stride 225000 \
                        --pool-sizes 40 \
                        --separator-char tab \
                        --out-dir ${OUTPUT_DIR} \
                        --file-prefix ${FILENAME%.sync}_250k_
done
