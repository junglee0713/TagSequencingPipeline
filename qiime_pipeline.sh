#!/bin/bash
set -x
set -e
set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 WORK_DIR MAPPING_FP LANE_NUM"
    exit 1
fi

WORK_DIR="$1"
MAPPING_FP="$2"
LANE_NUM="$3"

DATA_DIR="${WORK_DIR}/data_files"
JOINED_DIR="${WORK_DIR}/qiime_joined"
LIBRARY_DIR="${WORK_DIR}/library"

QIIME_PARAMS_FP="/home/tanesc/qiime_parameters.txt"

#MAPPING_FP="${WORK_DIR}/sampleSheet.tsv"
FWD="${DATA_DIR}/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_R1_001.fastq"
REV="${DATA_DIR}/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_R2_001.fastq"
IDX="${DATA_DIR}/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_I1_001.fastq"

PROJECTS=$(awk 'NR>1 {print $NF}' $MAPPING_FP | uniq)

#join_paired_ends.py \
#    -f $FWD \
#    -r $REV \
#    -b $IDX \
#    --min_overlap 35 \
#    --perc_max_diff 15 \
#    -o $JOINED_DIR

#split_libraries_fastq.py \
#    -i "${JOINED_DIR}/fastqjoin.join.fastq" \
#    -b "${JOINED_DIR}/fastqjoin.join_barcodes.fastq" \
#    -m $MAPPING_FP \
#    -q19 \
#    -o $LIBRARY_DIR \
#    --rev_comp_mapping_barcodes \
#    --barcode_type golay_12 \
#    --max_barcode_errors=2

## for each project pick otus, calculate beta diversities  
for PROJECT in $PROJECTS
do
    echo $PROJECT
    PROJECT_DIR="${WORK_DIR}/{PROJECT}"
    mkdir $PROJECT_DIR

    PROJECT_MAPPING_FP="${PROJECT_DIR}/${PROJECT}_sampleSheet.tsv"
    

    PROJECT_SAMPLE_ID="${PROJECT_DIR}/{PROJECT}_sampleNames.txt"
    $(cut -f 1 "$PROJECT_MAPPING_FP") > $PROJECT_SAMPLE_ID
    
    PROJECT_LIBRARY_DIR="${PROJECT_DIR}/library"
    PROJECT_OTU_DIR="${PROJECT_DIR}/otu"
    PROJEC_BETA_DIR="${PROJECT_DIR}/beta_diversity"
    
    #filter_fasta.py \
#	-f "{LIBRARY_DIR}/seqs.fna" \
#	-o "${PROJECT_LIBRARY_DIR}/seqs.fna" \
#	--sample_id_fp $PROJECT_SAMPLE_ID
    
#    cp "{LIBRARY_DIR}/split_library_log" $PROJECT_LIBRARY_DIR
    
#    pick_de_novo_otus.py \
#	--input_fp "${PROJECT_LIBRARY_DIR}/seqs.fna" \
#	--output_dir $PROJECT_OTU_DIR \
#	-p $QIIME_PARAMS_FP

#    biom convert \
#	-i "${PROJECT_OTU_DIR}/otu_table.biom"  \
#	-o "${PROJECT_OTU_DIR}/otu_table.txt" \
#	--to-tsv \
#	-m $PROJECT_MAPPING_FP \
#	--header-key=taxonomy \
#	--output-metadata-id="Consensus Lineage" \
#	--process-obs-metadata=taxonomy

#    beta_diversity_through_plots.py \
#	--otu_table_fp "${PROJECT_OTU_DIR}/otu_table.biom" \
#	--mapping_fp $MAPPING_FP \
#	--parameter_fp $QIIME_PARAMS_FP \
#	--output $PROJECT_BETA_DIR \
#	--tree_fp "${PROJECT_OTU_DIR}/rep_set.tre"
done
