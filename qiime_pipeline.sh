#!/bin/bash
set -x
set -e
set -u

if [ $# -ne 1 ]; then
    echo "Usage: $0 WORK_DIR"
    exit 1
fi

WORK_DIR="$1"
DATA_DIR="${WORK_DIR}/data_files"
JOINED_DIR="${WORK_DIR}/qiime_joined"
LIBRARY_DIR="${WORK_DIR}/library"

QIIME_PARAMS_FP="/home/tanesc/qiime_parameters.txt"

MAPPING_FP="${WORK_DIR}/sampleSheet.tsv"
FWD="${DATA_DIR}/Undetermined_S0_L001_R1_001.fastq.gz"
REV="${DATA_DIR}/Undetermined_S0_L001_R2_001.fastq.gz"
IDX="${DATA_DIR}/Undetermined_S0_L001_I1_001.fastq.gz"

PROJECTS=$(awk '{print $NF}' friedman_qiime_mapping.tsv | uniq)

join_paired_ends.py \
    -f $FWD \
    -r $REV \
    -b $IDX \
    --min_overlap 35 \
    --perc_max_diff 15 \
    -o $JOINED_DIR

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

    ########### where do i make the project specific sample sheet
    PROJECT_MAPPING_FP="" ######

    PROJECT_SAMPLE_ID="${PROJECT}_sampleNames.txt" ######### look into this!!!!!!!!!!
    $(cut -f 1 "$PROJECT_MAPPING_FP") > $PROJECT_SAMPLE_ID

    PROJECT_DIR="${WORK_DIR}/{PROJECT}"
    mkdir $PROJECT_DIR
    
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
