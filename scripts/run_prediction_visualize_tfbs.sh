#!/bin/bash

###############################################################################
# Script: run_prediction_visualize_tfbs.sh
# Description:
#   Runs DNABERT predictions with attention visualization on variant-wise
#   sequences for each TFBS model. This script reads TFBS model folder names
#   from a text file and processes each using its best fine-tuned model.
#
# Input:
#   - TFBS names listed in: $tfbs_file_path
#   - dev.tsv in: $base_dir_path/<TFBS>/variant_wise/
#   - Model checkpoints in: $base_model_path/<TFBS>/e-*/<checkpoint>
#
# Output:
#   - Prediction and attention scores saved in: Prediction_result/
#
# Author: Pratik Dutta
# Date: 04/30/2025
###############################################################################

# Set DNABERT configuration parameters
export KMER=6
export CANCER_TYPE="Brain"
export ARCHITECTURE="GDC_${CANCER_TYPE}_TFBS_Prediction_Model"
export CLASSES_NAME="NonTFBS_TFBS"

# Assign GPU devices (adjust as needed)
export CUDA_VISIBLE_DEVICES=5,6,7,0,1

# Define key base directories
base_dir_path="/home/pdutta/Data/Cancer_wiseGDC/New_data/Brain/Generated_files/DNABERT_Data/Somatic/300bp_TFBS_Agg"
base_model_path="/home/pdutta/Data/TFBS_Best_finetuned_models/Best_models"
tfbs_file_path="/home/pdutta/Github/Postdoc/DNABERT_data_processing/TFBS/top_tfbs_visualization_list/tfbs_visualize.txt"

# Move to root of repository (if needed)
cd ../..
cd ..

# Loop through each TFBS name from the input file
while IFS= read -r folder_name; do
    echo "------------------------------------------------------------"
    echo "📂 Processing TFBS model: $folder_name"
    
    # Path to TFBS input folder (containing dev.tsv files)
    folder="${base_dir_path}/${folder_name}/"
    
    # Locate the best model directory: e-*/checkpoint/
    model_e_dir=$(ls -d ${base_model_path}/${folder_name}/e-* 2>/dev/null | head -n 1)
    if [ -z "$model_e_dir" ]; then
        echo "❌ No 'e-*' model directory found for $folder_name. Skipping..."
        continue
    fi

    model_subdir=$(ls -d ${model_e_dir}/* 2>/dev/null | head -n 1)
    if [ -z "$model_subdir" ]; then
        echo "❌ No subdirectory inside $model_e_dir. Skipping..."
        continue
    fi

    export MODEL_PATH=$model_subdir
    echo "✅ Using model checkpoint: $MODEL_PATH"

    # Path to patient-wise variant input
    parent_dir_path="${folder}variant_wise"
    echo "🔍 Searching dev.tsvs in: $parent_dir_path"

    # Loop through each patient subfolder in variant_wise
    for subfolder in $(ls -d "${parent_dir_path}"); do
        echo "📄 Processing patient folder: $(basename "$subfolder")"

        # Define paths for prediction output
        export DATA_PATH=$subfolder
        export PREDICTION_PATH=$subfolder/Prediction_result
        export SUMMARY_PATH=$subfolder/Prediction_result
        export TB_PATH=$subfolder/Prediction_result

        # Run DNABERT prediction and generate attention scores
        python run_finetune_WANDB.py \
            --model_type dna \
            --tokenizer_name=dna$KMER \
            --model_name_or_path $MODEL_PATH \
            --architecture $ARCHITECTURE \
            --task_name dnaprom \
            --classes_name $CLASSES_NAME \
            --do_visualize \
            --data_dir $DATA_PATH  \
            --visualize_data_dir $PREDICTION_PATH \
            --visualize_models $KMER \
            --max_seq_length 300 \
            --per_gpu_pred_batch_size=100 \
            --output_dir $MODEL_PATH \
            --summary_dir $SUMMARY_PATH \
            --predict_dir $PREDICTION_PATH \
            --tb_log_dir $TB_PATH \
            --wandb_tags GDC_Genomic_VCF TFBS $CANCER_TYPE $(basename "$subfolder") $(basename "$folder") \
            --n_process 66
    done

done < "$tfbs_file_path"

echo "🎉 All TFBS model predictions completed!"
