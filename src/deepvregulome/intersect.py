import os
import pandas as pd
import numpy as np
import pysam
import pickle
import subprocess
import pybedtools
import shutil
import multiprocessing
import tempfile
import yaml

# Load configuration from config.yaml
with open("config.yaml") as f:
    config = yaml.safe_load(f)

# Extracting paths from configuration
vcf_folder_path = config["vcf_folder_path"]
chip_seq_base_path = config["chip_seq_base_path"]
output_base_path = config["output_base_path"]
reference_genome_path = config["reference_genome_path"]
vcf_error_log_dir = config["vcf_error_log_dir"]
pybedtools.helpers.set_tempdir(config["pybedtools_temp"])

# Load reference genome for sequence extraction (not yet used here)
reference_fasta = pysam.FastaFile(reference_genome_path)

# Delete all contents from temporary directory (e.g. pybedtools tempdir)
def clear_temp_directory(temp_dir):
    for filename in os.listdir(temp_dir):
        file_path = os.path.join(temp_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

# Recursively collect all VCF files, excluding logs directory
def get_vcf_gz_files_except_logs(root_folder):
    all_files = []
    for dirpath, dirnames, filenames in os.walk(root_folder):
        if 'logs' in dirnames:
            dirnames.remove('logs')
        for filename in filenames:
            if filename.endswith('.vcf.gz'):
                all_files.append(os.path.join(dirpath, filename))
    return all_files

# Reindex VCF using tabix for fast querying (optional utility)
def reindex_vcf(vcf_path):
    try:
        subprocess.run(['tabix', '-p', 'vcf', vcf_path], check=True)
        print(f"Re-indexed successfully: {vcf_path}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to re-index {vcf_path}: {e}")

# Convert a VCF file to a DataFrame with coordinates and metadata
def vcf_to_dataframe(vcf_path):
    try:
        vcf_file = pysam.VariantFile(vcf_path)
    except ValueError as e:
        print(f"Failed to open file {vcf_path}: {e}")
        os.makedirs(vcf_error_log_dir, exist_ok=True)
        with open(os.path.join(vcf_error_log_dir, vcf_path.split("/")[-2] + "_error_log.txt"), "a") as error_log:
            error_log.write(f"{vcf_path}\n")
        return None

    data = []
    columns = ["CHROM", "START_POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + list(vcf_file.header.info.keys())
    for record in vcf_file:
        basic_data = [record.chrom, record.pos, record.id, record.ref, ','.join(str(alt) for alt in record.alts), record.qual, record.filter.keys()[0] if record.filter.keys() else 'PASS']
        row_data = [record.info.get(key) for key in vcf_file.header.info.keys()]
        data.append(basic_data + row_data)

    df = pd.DataFrame(data, columns=columns)
    df["START_POS"] = df["START_POS"] - 1  # Convert to 0-based
    df.insert(2, 'END_POS', df["START_POS"] + df['REF'].str.len())
    vcf_file.close()
    return df

# Main function that intersects each TFBS model with all VCFs
def process_tfbs_data(sub_df, temp_dir, files):
    for _, row in sub_df.iterrows():
        tf_name = row[0]
        file_path = f"{chip_seq_base_path}/{tf_name}/300bp_unique_balanced.bed"
        print(file_path)
        tfbs_output_path = f"{output_base_path}/{tf_name}"
        os.makedirs(tfbs_output_path, exist_ok=True)

        df_tfbs = pd.read_csv(file_path, sep="\t")
        TFBS_bed = pybedtools.BedTool.from_dataframe(df_tfbs)
        intersected_vcf_data = {}
        df_statistics = pd.DataFrame(columns=["filename", "Patient_ID", 'work_flow', 'VCF_instance', 'VCF_feature', 'Intersected_instances', 'VCF_column_names'])

        for file_path in files:
            file_name = os.path.basename(file_path)
            parts = file_name.split('.')
            df_vcf = vcf_to_dataframe(file_path)
            if df_vcf is None:
                continue

            vcf_bed = pybedtools.BedTool.from_dataframe(df_vcf)
            intersect_temp_path = os.path.join(temp_dir, "intersect_temp.bed")
            TFBS_bed.intersect(vcf_bed, wa=True, wb=True, output=intersect_temp_path)

            column_list = df_tfbs.columns.to_list() + df_vcf.columns.to_list()
            df_intersection = pybedtools.BedTool(intersect_temp_path).to_dataframe(names=column_list)
            df_intersection = df_intersection[(df_intersection['REF'].str.len() < 10) & (df_intersection['ALT'].str.len() < 10)]

            patient_ID, work_flow = parts[0], parts[2]
            df_statistics.loc[len(df_statistics)] = [file_name, patient_ID, work_flow, df_vcf.shape[0], df_vcf.shape[1], df_intersection.shape[0], list(df_vcf.columns)]
            intersected_vcf_data[f"{patient_ID}_{work_flow}"] = df_intersection

        # Save intersection stats and intersected data per TFBS model
        df_statistics.to_csv(os.path.join(tfbs_output_path, "VCF_statistics.tsv"), sep="\t", index=False)
        with open(os.path.join(tfbs_output_path, "intersected_vcf_data.pkl"), "wb") as file:
            pickle.dump(intersected_vcf_data, file)

# Split a DataFrame into N parts (used for multiprocessing)
def split_dataframe(df, num_splits):
    return np.array_split(df, num_splits)

# Wrapper to call TFBS intersection for a batch in a temp dir
def parallel_process_tfbs_data(sub_df, files):
    with tempfile.TemporaryDirectory(dir=config["pybedtools_temp"]) as temp_dir:
        process_tfbs_data(sub_df, temp_dir, files)

# Entrypoint: given a TFBS model list directory, run the intersection pipeline
def run_tfbs_intersection_pipeline(tfbs_list_dir, start=4, end=46, num_processes=5):
    files = get_vcf_gz_files_except_logs(vcf_folder_path)
    print(f"Total VCFs: {len(files)}")

    for file_index in range(start, end):
        print(f"Processing TFBS list file: output_{file_index}.txt")
        tfbs_list_path = os.path.join(tfbs_list_dir, f"output_{file_index}.txt")
        tfbs_df = pd.read_csv(tfbs_list_path, sep="\t", header=None)

        sub_dfs = split_dataframe(tfbs_df, num_processes)
        pool = multiprocessing.Pool(processes=num_processes)
        pool.starmap(parallel_process_tfbs_data, [(sub_df, files) for sub_df in sub_dfs])
        pool.close()
        pool.join()

        clear_temp_directory(config["pybedtools_temp"])
        print(f"Cleaned up temp after processing file index {file_index}")