import os
import pandas as pd
import pysam
import gzip
import pickle, subprocess
import pybedtools
print(pybedtools.__file__)
pybedtools.helpers.set_tempdir('/home/pdutta/pybedtools_temp') ## Set a temporary folder


def get_vcf_gz_files_except_logs(root_folder):
    all_files = []

    # Walk through the directory tree
    for dirpath, dirnames, filenames in os.walk(root_folder):
        # If "logs" is in dirnames, remove it to avoid traversing it
        if 'logs' in dirnames:
            dirnames.remove('logs')

        # Add only the filenames with the extension .vcf.gz in the current directory to the all_files list
        for filename in filenames:
            if filename.endswith('.vcf.gz'):
                all_files.append(os.path.join(dirpath, filename))
    return all_files



def vcf_to_dataframe(file_info):
    """
    Convert a .vcf.gz file into a pandas DataFrame.

    Parameters:
    - file_info (tuple): contains (vcf_path, patient_id)

    Returns:
    - pd.DataFrame: VCF data as a DataFrame with patient_id column
    """
    vcf_path, patient_id = file_info
    print(patient_id)

    try:
        # Open the VCF file
        vcf_file = pysam.VariantFile(vcf_path)
    except ValueError as e:
        print(f"Failed to open file {vcf_path}: {e}")
        with open("/data/private/pdutta_new/GDC_Cancer_Wise/New_data/Brain/" + vcf_path.split("/")[-2] + "_error_log.txt", "a") as error_log:
            error_log.write(f"{vcf_path}\n")
        return None

    # Extracting the data and the columns
    data = []
    columns = ["CHROM", "START_POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + list(vcf_file.header.info.keys())
    for record in vcf_file:
        basic_data = [record.chrom, record.pos, record.id, record.ref,
                      ','.join(str(alt) for alt in record.alts), record.qual,
                      record.filter.keys()[0] if record.filter.keys() else 'PASS']
        row_data = [record.info.get(key) for key in vcf_file.header.info.keys()]
        data.append(basic_data + row_data)

    df = pd.DataFrame(data, columns=columns)
    df["START_POS"] = df["START_POS"] - 1
    end = df["START_POS"] + df['REF'].str.len()
    df.insert(2, 'END_POS', end)

    # Add the "patient_id" column (optional )
    df['patient_id'] = patient_id

    # Close the VCF file
    vcf_file.close()

    return df


def process_files_in_parallel(files):
    # Set up multiprocessing
    with Pool(cpu_count() - 1) as pool:
        results = pool.map(vcf_to_dataframe, files)
    
    # Filter out any None results and concatenate the DataFrames
    return pd.concat([df for df in results if df is not None], ignore_index=True)