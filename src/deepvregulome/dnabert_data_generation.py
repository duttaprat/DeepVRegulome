import os, glob
import pandas as pd
import pysam
import numpy as np
import pickle
from multiprocessing import Pool, cpu_count

# === Configurable paths ===
cancer_type = "Brain"
intersected_base_path = f"/data/private/pdutta_new/GDC_Cancer_Wise/New_data/{cancer_type}/Generated_files/Intersected_Data/Somatic/300bp_TFBS"
reference_genome_path = "/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"
output_path = f"/data/private/pdutta_new/GDC_Cancer_Wise/New_data/{cancer_type}/Generated_files/DNABERT_Data/Somatic/TFBS_300bp_new"

os.makedirs(output_path, exist_ok=True)

# === Load top models list ===
top_models_df = pd.read_csv("/home/campus.stonybrook.edu/pdutta/Github/Postdoc/DNABERT_data_processing/TFBS/300bp_TFBS_accuracy_Stat.tsv", sep="\t")
top_models_df = top_models_df[top_models_df['eval_acc'] >= 0.85].iloc[100:200]

# === Convert sequence to space-separated k-mers ===
def seq2kmer(seq, k=6):
    return " ".join([seq[i:i+k] for i in range(len(seq)-k+1)])

# === Apply all mutations to a given reference sequence ===
def apply_mutations_absolute(ref_seq, chip_start, starts, ends, refs, alts):
    mutations = sorted(zip(starts, ends, refs, alts), key=lambda x: x[0])
    seq_list = list(ref_seq)
    offset = chip_start
    for start, end, ref, alt in mutations:
        rel_pos = start - offset
        if ''.join(seq_list[rel_pos:rel_pos+len(ref)]) == ref:
            seq_list[rel_pos:rel_pos+len(ref)] = list(alt)
    return ''.join(seq_list)

# === Load genome reference once per worker ===
def initialize_worker():
    global reference_fasta
    reference_fasta = pysam.FastaFile(reference_genome_path)
    print("##### Reference genome loaded in worker #####")

# === Process a dataframe of variants for a single TFBS-patient combo ===
def get_sequences(df):
    global reference_fasta
    try:
        data = []
        for idx, row in df.iterrows():
            chrom, ref_start, ref_end = row[0], row[1], row[2]
            variant_start, variant_end = row['START_POS'], row['END_POS']
            ref_nucleotide, alt = row['REF'], row['ALT']
            ref_seq = reference_fasta.fetch(chrom, ref_start, ref_end)
            variant_pos_start = variant_start - ref_start
            variant_pos_end = variant_end - ref_start

            if len(ref_nucleotide) < len(alt):  # Insertion
                delete_size = len(alt) - len(ref_nucleotide)
                alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:len(ref_seq)-delete_size]
            elif len(ref_nucleotide) > len(alt):  # Deletion
                insert_size = len(ref_nucleotide) - len(alt)
                extra_bases = reference_fasta.fetch(chrom, ref_end, ref_end+insert_size)
                alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:] + extra_bases
            else:  # SNV
                alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:]

            data.append({
                'chr': chrom, 'Chip_Seq_start': ref_start, 'Chip_Seq_end': ref_end,
                'varinat_start': variant_start, 'variant_end': variant_end,
                'ref_neucleotide': ref_nucleotide, 'alternative_neucleotide': alt,
                'reference_seq': ref_seq, 'alt_seq': alt_seq
            })

        new_df = pd.DataFrame(data).drop_duplicates().reset_index(drop=True)

        merged_list = list(zip(new_df['reference_seq'], new_df['alt_seq']))
        merged_list = [item.upper() for tup in merged_list for item in tup]
        df_kmer = pd.DataFrame(list(map(seq2kmer, merged_list)), columns=['Sequence'])
        df_kmer['Label'] = np.random.choice([0, 1], size=len(df_kmer))

        grouped_df = new_df.groupby(['chr', 'Chip_Seq_start', 'Chip_Seq_end']).agg({
            'varinat_start': list, 'variant_end': list,
            'ref_neucleotide': list, 'alternative_neucleotide': list,
            'reference_seq': 'first', 'alt_seq': lambda x: list(set(x))
        }).reset_index()

        grouped_df['mutated_sequence'] = grouped_df.apply(lambda row: apply_mutations_absolute(
            row['reference_seq'], row['Chip_Seq_start'], row['varinat_start'],
            row['variant_end'], row['ref_neucleotide'], row['alternative_neucleotide']), axis=1)

        merged_regionwise = list(zip(grouped_df['reference_seq'], grouped_df['mutated_sequence']))
        merged_regionwise = [item.upper() for tup in merged_regionwise for item in tup]
        df_kmer_region = pd.DataFrame(list(map(seq2kmer, merged_regionwise)), columns=['Sequence'])
        df_kmer_region['Label'] = np.random.choice([0, 1], size=len(df_kmer_region))

        return new_df, df_kmer, grouped_df, df_kmer_region

    except Exception as e:
        print(f"Error in get_sequences: {e}")
        raise

# === Main processing logic for each TFBS tag ===
def process_row(row):
    global reference_fasta
    try:
        print(f"Processing {row['tags']}")
        intersected_data = f"{intersected_base_path}/{row['tags']}/intersected_vcf_data.pkl"

        with open(intersected_data, "rb") as file:
            loaded_dictionary = pickle.load(file)

        dnabert_raw_data = {}
        dnabert_raw_data_regionwise = {}

        for key, value in loaded_dictionary.items():
            try:
                print(f"  → {key}")
                new_df, df_kmer, grouped_df, df_kmer_region = get_sequences(value)
                dnabert_raw_data[key] = new_df
                dnabert_raw_data_regionwise[key] = grouped_df

                patient_folder = f"{output_path}/{row['tags']}/Patient_wise"
                os.makedirs(f"{patient_folder}/variant_wise/{key}", exist_ok=True)
                os.makedirs(f"{patient_folder}/region_wise/{key}", exist_ok=True)
                df_kmer.to_csv(f"{patient_folder}/variant_wise/{key}/dev.tsv", sep="\t", index=False)
                df_kmer_region.to_csv(f"{patient_folder}/region_wise/{key}/dev.tsv", sep="\t", index=False)

            except Exception as e:
                print(f"Error: {e} in {row['tags']} → {key}")
                return row['tags']

        with open(f"{output_path}/{row['tags']}/variantwise_raw_vcf_data.pkl", "wb") as file:
            pickle.dump(dnabert_raw_data, file)
        with open(f"{output_path}/{row['tags']}/regionwise_raw_vcf_data.pkl", "wb") as file:
            pickle.dump(dnabert_raw_data_regionwise, file)

    except Exception as e:
        print(f"Global error in {row['tags']}: {e}")
        return row['tags']

    return None

# === Run multiprocessing ===
if __name__ == '__main__':
    print(f"Using {cpu_count()} CPUs for processing")
    with Pool(cpu_count() - 10, initializer=initialize_worker) as pool:
        missing_files = pool.map(process_row, [row for _, row in top_models_df.iterrows()])

    missing_files = [tag for tag in missing_files if tag is not None]
    with open(f"{output_path}/missing_files.txt", "w") as f:
        for tag in missing_files:
            f.write(f"{tag}\n")

    print(f"Missing files recorded at {output_path}/missing_files.txt")
