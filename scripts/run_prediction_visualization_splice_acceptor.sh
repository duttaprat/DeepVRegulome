#!/bin/bash

###############################################################################
# Script: run_prediction_visualize_splice_acceptor.sh
# Description:
#   Runs DNABERT predictions on somatic variant-region sequences associated
#   with Acceptor splice sites in GDC Brain Cancer samples. Each patient's
#   dev.tsv file (region-wise) is used as input. Model outputs and attention
#   scores are saved to Prediction_result folder inside each subfolder.
#
# Requires:
#   - dev.tsv files in: DNABERT_Data/Somatic/acceptor/Patient_wise/variant_wise/
#   - Trained model at: $MODEL_PATH
#   - DNABERT visualization script: run_finetune_WANDB.py
#
# Author: [Pratik Dutta]
# Date: [04/30/2025]
###############################################################################

export LC_ALL=en_US.utf8
export LANG=en_US.utf8
export KMER=6
export CANCER_TYPE="Brain"
export ARCHITECTURE="GDC_New_${CANCER_TYPE}_Acceptor_Prediction_Model"
export CLASSES_NAME="AcceptorSpliceSite_NonAcceptorSpliceSite"
export MODEL_PATH="/data/private/pdutta_new/DNABERT_best_models/Acceptor/logstep=50_bs=640_lr=0.0003_wp=0.1_dp=0.1_wd=0.0005_len=80_epoch=20.0"

# Set path to variant-wise dev.tsv files
parent_dir_path="/data/projects/GDC_Cancer_Wise/New_data/${CANCER_TYPE}/Generated_files/DNABERT_Data/Somatic/acceptor/Patient_wise/variant_wise"
echo "🔍 Scanning variant-wise patient folders in: $parent_dir_path"

# Loop through each patient subfolder
for subfolder in "${parent_dir_path}"/*; do
    echo "------------------------------------------------------------"
    echo "📂 Processing: $(basename "$subfolder")"
    
    export DATA_PATH="$subfolder"
    export PREDICTION_PATH="${subfolder}/Prediction_result"
    export SUMMARY_PATH="$PREDICTION_PATH"
    export TB_PATH="$PREDICTION_PATH"

    mkdir -p "$PREDICTION_PATH"

    python run_finetune_WANDB.py \
        --model_type dna \
        --tokenizer_name=dna${KMER} \
        --model_name_or_path "$MODEL_PATH" \
        --architecture "$ARCHITECTURE" \
        --task_name dnaprom \
        --classes_name "$CLASSES_NAME" \
        --do_visualize \
        --data_dir "$DATA_PATH" \
        --max_seq_length 80 \
        --per_gpu_pred_batch_size=512 \
        --output_dir "$MODEL_PATH" \
        --summary_dir "$SUMMARY_PATH" \
        --predict_dir "$PREDICTION_PATH" \
        --tb_log_dir "$TB_PATH" \
        --wandb_tags GDC_VCF_Somatic Acceptor_regionwise "$CANCER_TYPE" "$(basename "$subfolder")" \
        --n_process 36

    echo "✅ Done with: $(basename "$subfolder")"
done

echo "🎉 All predictions complete!"
