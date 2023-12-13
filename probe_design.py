import csv
import os
import re
import requests

from math import log2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns



def subset_domain_fasta(input_fasta, output_fasta, sequence_names):
    """
    Extracts sequences from an input FASTA file based on provided sequence names and writes them to an output FASTA file.

    Parameters:
    - input_fasta (str): Path to the input FASTA file.
    - output_fasta (str): Path to the output FASTA file where the extracted sequences will be saved.
    - sequence_names (list): List of sequence names to be extracted from the input FASTA file.

    Reads the input FASTA file, extracts sequences whose names match those provided in `sequence_names`, and writes
    these sequences to the output FASTA file. It also generates a warning message if any provided sequence names
    are not found in the input FASTA file.

    Examples:
    ```
    subset_domain_fasta("input.fasta", "output.fasta", ["seq1", "seq2"])
    ```
    """
    sequences = []
    unique_sequence_names = set(sequence_names)
    missing_sequences = list(unique_sequence_names)
    with open(input_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in unique_sequence_names:
                sequences.append(record)
                unique_sequence_names.remove(record.id)
                missing_sequences.remove(record.id)
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    if missing_sequences:
        print(
            f"WARNING: Sequence(s) not found in the input FASTA file: {', '.join(missing_sequences)}; \nFASTA written to {output_fasta}"
        )
    else:
        print(f"FASTA written to {output_fasta}")


def max_homopolymer_length(seq):
    max_length = 0
    current_length = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_length += 1
        else:
            current_length = 1

        max_length = max(max_length, current_length)

    return max_length


def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    gc_count = sequence.count("G") + sequence.count("C")
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    return gc_content


def calculate_entropy(sequence):
    """Calculate the entropy of a DNA sequence."""
    probabilities = [float(sequence.count(base)) / len(sequence) for base in "ACGT"]
    entropy = -sum(p * log2(p) if p > 0 else 0 for p in probabilities)
    entropy = entropy / log2(4)
    return entropy


def score_sequence_repetition(sequence, min_n=2, max_n=5):
    """This function scores the degree of repetitive sequences in an input DNA seq. Default nmer range is 2-5."""
    nmer_counts = {n: {} for n in range(min_n, max_n + 1)}

    for n in range(min_n, max_n + 1):
        for i in range(len(sequence) - n + 1):
            nmer = sequence[i : i + n]
            if nmer in nmer_counts[n]:
                nmer_counts[n][nmer] += 1
            else:
                nmer_counts[n][nmer] = 1

    repetitiveness_score = sum(
        sum(len(nmer) * count for nmer, count in counts.items() if count > 1)
        for counts in nmer_counts.values()
    )

    return repetitiveness_score


def nominate_probes(
    fasta_path,
    gc_diff_threshold=14,
    homopolymer_length_threshold=4,
    entropy_threshold=0.8,
    repetitive_sequence_score_threshold=125,
):
    """
    This function identifies and nominates probes from a given FASTA file based on specified criteria.

    Parameters:
    -----------
    fasta_path : str
        The path to the input FASTA file containing DNA sequences for probe design.

    gc_diff_threshold : int
        The threshold for the difference in GC content from average (58%). The default value is 14%.

    homopolymer_length_threshold : int
        The threshold for the maximum allowed homopolymer length in probe sequences. The default is 4 bases.

    entropy_threshold : int
        Threshold for the minimum sequence entropy as a surrogate for sequence diversity. See calculate_entropy() for more information. Default threshold is set to 0.8.

    repetitive_sequence_score_threshold : int
         Threshold for the maximum sequence repetitiveness calculated using score_sequence_repetition(). Numerates the number of 2-5 nmer repeats in LHS and RHS probes and aggregates to an int score.
         See score_sequence_repetition() for details. Default threshold is set to 125.

    Returns:
    --------
    Creates a directory probe_nomination in the directory of the FASTA input. Creates an annotation CSV with the measured probe 5' and 3' sites as well as the maximum homopolymer length and GC content.
    Also writes the probes that pass the threshold in FASTA format.

    Prints output file paths and the number of probe sites that passed the threshold.
    """
    fasta_handle = open(fasta_path)
    pass_sequences_records = []
    sequence_bins_count = {}
    sequences_with_passing_bins = 0

    base_filename = os.path.basename(fasta_path).rsplit(".", 1)[0]

    output_directory = os.path.join(os.getcwd(), "probe_nomination")
    os.makedirs(output_directory, exist_ok=True)

    for record in SeqIO.parse(fasta_handle, "fasta"):
        reverse_complement_sequence = record.seq.reverse_complement().upper()
        bin_count = 0

        for i in range(0, len(reverse_complement_sequence) - 50 + 1):
            bin_sequence = reverse_complement_sequence[i : i + 50]

            gc_content_lhs = round(calculate_gc_content(bin_sequence[:25]),3)
            gc_content_rhs = round(calculate_gc_content(bin_sequence[-25:]),3)
            lhs_max_homopolymer = max_homopolymer_length(bin_sequence[:25])
            rhs_max_homopolymer = max_homopolymer_length(bin_sequence[-25:])
            lhs_entropy = round(calculate_entropy(bin_sequence[:25]),3)
            rhs_entropy = round(calculate_entropy(bin_sequence[-25:]),3)
            lhs_repetitiveness = score_sequence_repetition(bin_sequence[:25])
            rhs_repetitiveness = score_sequence_repetition(bin_sequence[-25:])

            if (
                bin_sequence[24] == "T"
                and abs(gc_content_lhs - 58) <= gc_diff_threshold
                and abs(gc_content_rhs - 58) <= gc_diff_threshold
                and lhs_max_homopolymer <= homopolymer_length_threshold
                and rhs_max_homopolymer <= homopolymer_length_threshold
                and lhs_entropy >= entropy_threshold
                and rhs_entropy >= entropy_threshold
                and lhs_repetitiveness <= repetitive_sequence_score_threshold
                and rhs_repetitiveness <= repetitive_sequence_score_threshold
            ):
                seq_record = SeqRecord(
                    bin_sequence,
                    id=f"{record.id}_bin_{i + 1}",
                    description=f"50-bp bin starting at position {i + 1}",
                    annotations={
                        "lhs_gc_content": gc_content_lhs,
                        "rhs_gc_content": gc_content_rhs,
                        "lhs_max_homopolymer": lhs_max_homopolymer,
                        "rhs_max_homopolymer": rhs_max_homopolymer,
                        "lhs_entropy": lhs_entropy,
                        "rhs_entropy": rhs_entropy,
                        "lhs_repetitiveness": lhs_repetitiveness,
                        "rhs_repetitiveness": rhs_repetitiveness,
                        "lhs_sequence": str(bin_sequence[:25]),  
                        "rhs_sequence": str(bin_sequence[-25:]), 
                    },
                )
                pass_sequences_records.append(seq_record)
                bin_count += 1

        sequence_bins_count[record.id] = bin_count

        if bin_count > 0:
            sequences_with_passing_bins += 1

    output_directory = os.path.join(os.getcwd(), "probe_nomination")
    os.makedirs(output_directory, exist_ok=True)

    annotations_csv = os.path.join(output_directory, f"{base_filename}_annotations.csv")
    with open(annotations_csv, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            [
                "id",
                "description",
                "start",
                "end",
                "lhs_gc_content",
                "rhs_gc_content",
                "lhs_max_homopolymer",
                "rhs_max_homopolymer",
                "lhs_entropy",
                "rhs_entropy",
                "lhs_repetitiveness",
                "rhs_repetitiveness",
                "lhs_sequence",  # Include the columns for left-hand and right-hand sequences
                "rhs_sequence",
            ]
        )
        for seq_record in pass_sequences_records:
            start_position = int(seq_record.id.split("_")[-1])
            end_position = start_position + len(seq_record.seq) - 1
            csv_writer.writerow(
                [
                    seq_record.id,
                    seq_record.description,
                    start_position,
                    end_position,
                    seq_record.annotations.get("lhs_gc_content", ""),
                    seq_record.annotations.get("rhs_gc_content", ""),
                    seq_record.annotations.get("lhs_max_homopolymer", ""),
                    seq_record.annotations.get("rhs_max_homopolymer", ""),
                    seq_record.annotations.get("lhs_entropy", ""),
                    seq_record.annotations.get("rhs_entropy", ""),
                    seq_record.annotations.get("lhs_repetitiveness", ""),
                    seq_record.annotations.get("rhs_repetitiveness", ""),
                    seq_record.annotations.get("lhs_sequence", ""),  
                    seq_record.annotations.get("rhs_sequence", ""),
                ]
            )

    print(f"Annotations written to: {annotations_csv}")

    output_fasta = os.path.join(output_directory, f"{base_filename}_probes.fasta")
    SeqIO.write(pass_sequences_records, output_fasta, "fasta")
    print(f"Probes written to: {output_fasta}")

    fasta_handle.close()
    for record_id, bin_count in sequence_bins_count.items():
        print(
            f"Sequence {record_id} has {bin_count} potential probe binding sites passing threshold."
        )


def lhs_rhs_probe_binding_site_split(input_fasta, output_fasta):
    """This function ouputs the LHS and RHS probes as a FASTA from input probe-binding site nomination FASTA"""
    with open(input_fasta, "r") as input_file, open(output_fasta, "w") as output_file:
        for record in SeqIO.parse(input_file, "fasta"):
            # Extract first 25 bases (LHS) and last 25 bases (RHS)
            lhs_sequence = record.seq[:25]
            rhs_sequence = record.seq[-25:]

            # Create a new record with labeled sequences
            new_record = (
                f">{record.id}_LHS\n{lhs_sequence}\n>{record.id}_RHS\n{rhs_sequence}\n"
            )

            # Write the new record to the output file
            output_file.write(new_record)
    print(f"LHS and RHS probes written to: {output_fasta}")


def read_blast_results(file_path):
    """
    Reads a TSV file into a DataFrame, filters data, calculates columns,
    performs pattern matching, and extracts information.

    Parameters:
    - file_path (str): The path to the TSV file.
    - column_names (list): A list of column names for the DataFrame.

    Returns:
    - lhs_df, rhs_df (tuple of DataFrames): Two DataFrames containing processed data.
    """
    column_names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "stitle",
    ]
    blast_out = pd.read_csv(file_path, sep="\t", header=None, names=column_names)
    blast_out["off_target_binding"] = blast_out["length"] - blast_out["mismatch"]
    probe_max_offtarg = blast_out.groupby("qseqid")["off_target_binding"].idxmax()
    blast_out = blast_out.loc[probe_max_offtarg].copy()  

    pattern = r"\((.*?)\),"
    lhs_df = blast_out[
        blast_out["qseqid"].str.contains("LHS")
    ].copy()  
    rhs_df = blast_out[
        blast_out["qseqid"].str.contains("RHS")
    ].copy()  
    lhs_df.loc[:, "off_target_gene"] = lhs_df["stitle"].str.extract(
        pattern, expand=False
    )
    rhs_df.loc[:, "off_target_gene"] = rhs_df["stitle"].str.extract(
        pattern, expand=False
    )
    return lhs_df, rhs_df


def process_lhs_rhs(lhs_df, rhs_df, annotation_csv_path, output_csv_path):
    """
    Merge information from 'lhs_df' and 'rhs_df' DataFrames with an annotation table
    specified in a CSV file and save the merged DataFrame to a CSV file.

    Parameters:
    -----------
    lhs_df : pandas DataFrame
        DataFrame containing left-hand side data with a column 'qseqid' representing identifiers.
    rhs_df : pandas DataFrame
        DataFrame containing right-hand side data with a column 'qseqid' representing identifiers.
    annotation_csv_path : str
        Filepath to the CSV file containing annotation information with a column 'id'.
    output_csv_path : str
        Filepath to save the resulting merged DataFrame as a CSV file.

    Returns:
    --------
    None

    Explanation:
    ------------
    1. Reads the annotation information from 'annotation_csv_path' into a DataFrame named 'target_table'.
    2. Processes 'lhs_df' and 'rhs_df':
        - Extracts 'id' from 'qseqid' columns in 'lhs_df' and 'rhs_df' based on specific patterns
          ('(.+)_LHS' for 'lhs_df' and '(.+)_RHS' for 'rhs_df').
        - Merges 'target_table' with subset columns of 'lhs_df' and 'rhs_df' based on 'id'
          using a left join.
        - Renames the merged columns to indicate their relation to left-hand side ('lhs') and right-hand side ('rhs')
          information.
    3. Concatenates the resulting merged DataFrames ('lhs_merged' and 'rhs_merged') horizontally
       and saves the merged DataFrame to 'output_csv_path' as a CSV file, excluding the index column.
    """
    target_table = pd.read_csv(annotation_csv_path)

    lhs_df["id"] = lhs_df["qseqid"].str.extract(r"(.+)_LHS")
    lhs_merged = target_table.merge(
        lhs_df[["id", "length", "mismatch", "off_target_binding", "off_target_gene"]],
        on="id",
        how="left",
    ).rename(
        columns={
            "length": "lhs_off_target_length",
            "mismatch": "lhs_mismatch",
            "off_target_binding": "lhs_max_off_target_binding",
            "off_target_gene": "lhs_off_target_gene",
        }
    )

    rhs_df["id"] = rhs_df["qseqid"].str.extract(r"(.+)_RHS")
    rhs_merged = lhs_merged.merge(
        rhs_df[["id", "length", "mismatch", "off_target_binding", "off_target_gene"]],
        on="id",
        how="left",
    ).rename(
        columns={
            "length": "rhs_off_target_length",
            "mismatch": "rhs_mismatch",
            "off_target_binding": "rhs_max_off_target_binding",
            "off_target_gene": "rhs_off_target_gene",
        }
    )

    rhs_merged.to_csv(output_csv_path, index=False)


def calculate_priority_score(anno):
    """
    Calculate priority scores based on probe GC content, maximum off-target binding, and sequence bias (entropy).

    Parameters:
    - anno (DataFrame): Annotation DataFrame which is the output of process_lhs_rhs().

    Returns:
    - anno (DataFrame): DataFrame with added 'priority_score' column.
    """
    anno["lhs_max_off_target_binding"].fillna(0, inplace=True)
    anno["rhs_max_off_target_binding"].fillna(0, inplace=True)

    anno["priority_score"] = 200 / (
        ((abs(anno["lhs_gc_content"] - 50) + abs(anno["rhs_gc_content"] - 50)) / 2)
        * (
            1
            + (anno["lhs_max_off_target_binding"] / 10) ** 2
            + (anno["rhs_max_off_target_binding"] / 10) ** 2
        )
        / anno[["lhs_entropy", "rhs_entropy"]].min(axis=1)
    )

    anno["target"] = anno["id"].str.split("_").str[0]

    return anno


def nominate_non_overlapping_sets(group):
    """
    Generate sets of non-overlapping indices for probe pairs in an anno object. Additionally, excludes probe pairs that have the same off-target gene.

    Parameters:
    - group (DataFrameGroupBy): DataFrame grouped by a specific criterion, in this context the probe target.

    Returns:
    - sets (list): List containing sets of non-overlapping indices for each target group.
    """
    sets = []
    for idx, row in group.iterrows():
        non_overlapping = [idx]

        for _, next_row in group.iterrows():
            if len(non_overlapping) >= 3:
                break
            if idx != next_row.name:
                overlaps = False
                for nominated_index in non_overlapping:
                    if (
                        group.loc[nominated_index, "start"] - next_row["start"] <= 50
                    ) and (next_row["end"] - group.loc[nominated_index, "end"] <= 50):
                        overlaps = True
                        break
                if (
                    not overlaps
                    and (row["lhs_off_target_gene"] != row["rhs_off_target_gene"])
                    and (
                        next_row["lhs_off_target_gene"]
                        != next_row["rhs_off_target_gene"]
                    )
                ):
                    non_overlapping.append(next_row.name)

        sets.append(non_overlapping)

    return sets


def select_highest_priority_sets(sets, group):
    """
    Select sets with the highest sum of priority scores.

    Parameters:
    - sets (list): List of sets of indices generated on a nominated probe list from nominate_non_overlapping_sets().
    - group (DataFrame): DataFrame containing the annotation data.

    Returns:
    - nominated_set (list): Set of indices with the highest sum of priority scores.
    """
    max_sum_3 = -np.inf
    max_sum_2 = -np.inf
    max_sum_1 = -np.inf

    nominated_set_3 = []
    nominated_set_2 = []
    nominated_set_1 = []

    for combination in sets:
        total_score_indices = []
        if isinstance(combination, list):
            for subset in combination:
                if isinstance(subset, int):
                    total_score_indices.append(subset)
                else:
                    total_score_indices.extend(subset)
        else:
            total_score_indices.append(combination)

        total_score_df = group.loc[total_score_indices]

        if len(total_score_df) == 3:
            total_score = total_score_df["priority_score"].sum()
            if total_score > max_sum_3:
                max_sum_3 = total_score
                nominated_set_3 = total_score_indices

        elif len(total_score_df) == 2 and len(total_score_indices) > 1:
            total_score = total_score_df["priority_score"].sum()
            if total_score > max_sum_2:
                max_sum_2 = total_score
                nominated_set_2 = total_score_indices

        elif len(total_score_df) == 1 and len(total_score_indices) > 0:
            total_score = total_score_df["priority_score"].sum()
            if total_score > max_sum_1:
                max_sum_1 = total_score
                nominated_set_1 = total_score_indices

    if nominated_set_3:
        return nominated_set_3
    elif nominated_set_2:
        return nominated_set_2
    elif nominated_set_1:
        return nominated_set_1
    else:
        return None


def nominate_top_probe_set(input_csv, output_csv):
    """
    Process an annotation CSV file, calculate priority scores, select non-overlapping sets for each target,
    and store the highest priority set into a new CSV file.
    Expects the input CSV file to have the following columns: 
    'id', 'start', 'end', 'lhs_off_target_gene', 'rhs_off_target_gene', 'lhs_gc_content','rhs_gc_content', 'lhs_entropy', and 'rhs_entropy'.

    Parameters:
    - input_csv (str): Path to the input annotation CSV file.
    - output_csv (str): Path to the output CSV file to store the selected probe sets.

    Returns:
    - CSV file containing the highest priorirty probe sets.
    """
    annotation_df = pd.read_csv(input_csv)

    annotated_data = calculate_priority_score(annotation_df)
    annotated_data["target"] = annotated_data["id"].str.split("_").str[0]

    grouped_data = annotated_data.groupby("target")

    selected_sets = []

    for target, group in grouped_data:
        non_overlap_sets = nominate_non_overlapping_sets(group)

        highest_priority_set = select_highest_priority_sets(non_overlap_sets, group)

        selected_sets.append(highest_priority_set)
        if len(highest_priority_set) == 0:
            print(f"ERROR: No probes nominated for target {target}")
        elif len(highest_priority_set) == 1 or len(highest_priority_set) == 2:
            print(
                f"WARNING: Only {len(highest_priority_set)} probe(s) nominated for target {target}"
            )

    selected_sets_flat = [
        idx for sublist in selected_sets if sublist for idx in sublist
    ]
    annotated_data.loc[selected_sets_flat].to_csv(output_csv, index=False)


def construct_custom_flex_probe_opool(nominated_probe_csv, pool_name, lhs_r2_adaptor, n_barcodes, directory = 'opools'):
    """
    Constructs an Excel file containing probe sequences and pool names from a nominated probe CSV file. This xlsx should be able to be uploaded to IDT for ordering.

    Parameters:
    nominated_probe_csv (str): Path to the nominated probe CSV file from nominate_top_probe_set().
    output_order_name (str): Path to the output Excel file to be generated.
    pool_name_prefix (str): Prefix to be used in the pool names.
    n_barcodes (int): Number of barcode sequences to consider. Typically 4 or 16. 

    Returns:
    None
    """
    data = pd.read_csv(nominated_probe_csv)
    sequences = []
    barcode_ids = [
        'BC001', 'BC002', 'BC003', 'BC004', 'BC005', 'BC006', 'BC007', 'BC008',
        'BC009', 'BC010', 'BC011', 'BC012', 'BC013', 'BC014', 'BC015', 'BC016'
    ]
    barcode_sequences = [
        'ACTTTAGG', 'AACGGGAA', 'AGTAGGCT', 'ATGTTGAC', 'ACAGACCT', 'ATCCCAAC', 'AAGTAGAG', 'AGCTGTGA',
        'ACAGTCTG', 'AGTGAGTG', 'AGAGGCAA', 'ACTACTCA', 'ATACGTCA', 'ATCATGTG', 'AACGCCGA', 'ATTCGGTT'
    ]
    
    if n_barcodes:
        n_barcodes = min(n_barcodes, len(barcode_sequences))
    else:
        n_barcodes = len(barcode_sequences)
    
    for index, row in data.iterrows():
        lhs_sequence = row['lhs_sequence']
        rhs_sequence = row['rhs_sequence']
        
        lhs_seq =  lhs_r2_adaptor + lhs_sequence
        lhs_pool = pool_name + "_lhs"
        sequences.append({'Pool name': lhs_pool, 'Sequence': lhs_seq})

        for i in range(n_barcodes):
            seq = "/5Phos/" + rhs_sequence + "ACGCGGTTAGCACGTANN" + barcode_sequences[i] + "CGGTCCTAGCAA"
            pool = barcode_ids[i] + "_" + pool_name + "_rhs" 
            sequences.append({'Pool name': pool, 'Sequence': seq})
    
    oligo_pool = pd.DataFrame(sequences)
    oligo_pool = oligo_pool.sort_values(by='Pool name')
    
    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        oligo_pool.to_excel(writer, index=False)
    
 

def gene_to_ensembl_wta_filter(gene_names, input_csv = 'references/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'):
    """
    Fetches Ensembl IDs for a list of gene names and filters data from a CSV file based on matching Ensembl IDs.

    Args:
    - gene_names (list): List of gene names to fetch Ensembl IDs for.
    - input_csv (str, optional): Path to the CSV file containing gene data.
      Defaults to 'reference/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'.

    Returns:
    - filtered_data (pandas DataFrame): Filtered data from the CSV based on matched Ensembl IDs.
      If some gene names are not found in the Ensembl database, prints a message and returns filtered data for found genes.

    Note:
    - Uses Ensembl REST API to fetch Ensembl IDs for given gene names.
    - Retrieves data from the CSV file and filters it based on Ensembl IDs.
    - Prints messages for any missing Ensembl IDs or failed data retrieval from the Ensembl database.
    """
    ensembl_ids = {}
    base_url = "https://rest.ensembl.org/lookup/symbol/human/"

    for gene_name in gene_names:
        endpoint = f"{base_url}{gene_name}?"
        response = requests.get(endpoint, headers={"Content-Type": "application/json"})

        if response.ok:
            data = response.json()
            if 'id' in data:
                ensembl_ids[gene_name] = data['id']
            else:
                print(f"Ensembl ID not found for gene name: {gene_name}")
        else:
            print(f"Failed to retrieve data for gene name: {gene_name}. Check your input or try again later.")

    data = pd.read_csv(input_csv, names=['gene_id', 'probe_seq', 'probe_id', 'included', 'region'], skiprows=6)

    missing_genes = [gene for gene in gene_names if gene not in ensembl_ids]
    
    if missing_genes:
        print(f"Gene names not found in CSV: {', '.join(missing_genes)}")
        filtered_data = data[data['gene_id'].isin(ensembl_ids.values())]
        return filtered_data  

    filtered_data = data[data['gene_id'].isin(ensembl_ids.values())]
    return filtered_data


def generate_flex_sgrna_opool(a_position_protospacers, b_position_protospacers, pool_name, n_barcodes, directory = 'opools'):
    """
    From a list of A and B protospacers generates flex probes.

    Args:
    - a_position_protospacers (str): Sequence for position A protospacers.
    - b_position_protospacers (str): Sequence for position B protospacers.
    - pool_name (str): Name of the pool.
    - n_barcodes (int): Number of barcode sequences.

    Output:
    - Creates an Excel file named YYMMDD_pool_name.xlsx containing probe orders in directory opools.
    """
    sgrna_a_rhs_constant_region = 'CTGAAAC'
    sgrna_b_rhs_constant_region = 'AAAC'
    partial_capture_seqeuence_1 = 'CGGTCCTAGCAA'

    parital_truseq_r2 = "CAGACGTGTGCTCTTCCGATCT"
    lhs_constant_sequence = "ACGCGGTTAGCACGTANN"
    sgrna_a_lhs_constant_region = "CTTGCTATGCACTCTTGTGCTTAGCT"
    sgrna_b_lhs_constant_region = "GCTATGCTGTTTCCAGCTTAGCTCTT"

    a_seq = a_position_protospacers.apply(Seq)
    b_seq = b_position_protospacers.apply(Seq)

    pool_names_rhs = [f"RHS_{pool_name}"] * (len(a_seq) * 2)
    sequences_rhs = []

    for a, b in zip(a_seq, b_seq):
        sequences_rhs.extend([
            f"/5Phos/{sgrna_a_rhs_constant_region}{a.reverse_complement()}{partial_capture_seqeuence_1}",
            f"/5Phos/{sgrna_b_rhs_constant_region}{b.reverse_complement()}{partial_capture_seqeuence_1}"
        ])

    rhs_probe_order = pd.DataFrame({
        'Pool name': pool_names_rhs,
        'Sequence': sequences_rhs
    })

    pool_names_lhs = []
    sequences_lhs = []
    barcode_ids = [
        'BC001', 'BC002', 'BC003', 'BC004', 'BC005', 'BC006', 'BC007', 'BC008',
        'BC009', 'BC010', 'BC011', 'BC012', 'BC013', 'BC014', 'BC015', 'BC016'
    ]
    barcode_sequences = [
        'ACTTTAGG', 'AACGGGAA', 'AGTAGGCT', 'ATGTTGAC', 'ACAGACCT', 'ATCCCAAC', 'AAGTAGAG', 'AGCTGTGA',
        'ACAGTCTG', 'AGTGAGTG', 'AGAGGCAA', 'ACTACTCA', 'ATACGTCA', 'ATCATGTG', 'AACGCCGA', 'ATTCGGTT'
    ]

    for i in range(n_barcodes):
        pool_names_lhs.extend([
            f"{barcode_ids[i]}_LHS_{pool_name}",
            f"{barcode_ids[i]}_LHS_{pool_name}"
        ])

        sequences_lhs.extend([
            f"{parital_truseq_r2}{lhs_constant_sequence}{barcode_sequences[i]}{sgrna_a_lhs_constant_region}",
            f"{parital_truseq_r2}{lhs_constant_sequence}{barcode_sequences[i]}{sgrna_b_lhs_constant_region}"
        ])

    lhs_probe_order = pd.DataFrame({
        'Pool name': pool_names_lhs,
        'Sequence': sequences_lhs
    })

    oligo_pool = pd.concat([rhs_probe_order, lhs_probe_order])

    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        oligo_pool.to_excel(writer, index=False)
        

def generate_flex_gene_opool(gene_list, pool_name, n_barcodes, directory = 'opools'):
    """
    From a list of A and B protospacers generates flex probes.

    Args:
    - a_position_protospacers (str): Sequence for position A protospacers.
    - b_position_protospacers (str): Sequence for position B protospacers.
    - pool_name (str): Name of the pool.
    - n_barcodes (int): Number of barcode sequences.

    Output:
    - Creates an Excel file named YYMMDD_pool_name.xlsx containing probe orders in directory opools.
    """
  
    partial_capture_sequence_1 = 'CGGTCCTAGCAA'
    rhs_consatant_sequence = 'ACGCGGTTAGCACGTANN'
    parital_truseq_r2 = "CAGACGTGTGCTCTTCCGATCT"

    sequences = []
    
    barcode_ids = [
        'BC001', 'BC002', 'BC003', 'BC004', 'BC005', 'BC006', 'BC007', 'BC008',
        'BC009', 'BC010', 'BC011', 'BC012', 'BC013', 'BC014', 'BC015', 'BC016'
    ]
    barcode_sequences = [
        'ACTTTAGG', 'AACGGGAA', 'AGTAGGCT', 'ATGTTGAC', 'ACAGACCT', 'ATCCCAAC', 'AAGTAGAG', 'AGCTGTGA',
        'ACAGTCTG', 'AGTGAGTG', 'AGAGGCAA', 'ACTACTCA', 'ATACGTCA', 'ATCATGTG', 'AACGCCGA', 'ATTCGGTT'
    ]
    
    for index, row in gene_list.iterrows():
        lhs_sequence = row['probe_seq'][:25]
        rhs_sequence = row['probe_seq'][-25:]
        
        lhs_seq =  parital_truseq_r2 + lhs_sequence
        lhs_pool = pool_name + "_lhs"
        sequences.append({'Pool name': lhs_pool, 'Sequence': lhs_seq})

        for i in range(n_barcodes):
            seq = "/5Phos/" + rhs_sequence + rhs_consatant_sequence + barcode_sequences[i] + partial_capture_sequence_1
            pool = barcode_ids[i] + "_" + pool_name + "_rhs" 
            sequences.append({'Pool name': pool, 'Sequence': seq})

    sequences_df = pd.DataFrame(sequences)
    sequences_df = sequences_df.sort_values(by='Pool name')
    
    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        sequences_df.to_excel(writer, index=False)
        
        
        
def score_plot_probes(lhs, rhs, title):
    """
    Calculate metrics and plot histograms for the given left-hand side (LHS) and right-hand side (RHS) sequences.

    Parameters:
    - lhs (Pandas Series): Pandas Series containing left-hand side sequences.
    - rhs (Pandas Series): Pandas Series containing right-hand side sequences.
    - title (str): Title for the plotted histograms.

    This function calculates metrics including Max Homopolymer Length, GC Content, Entropy,
    and Repetitiveness Scores for both LHS and RHS sequences. It then displays histograms
    for each metric comparison between LHS and RHS sequences using Seaborn.

    The histograms include:
    - Distribution of Max Homopolymer Length
    - Distribution of GC Content
    - Distribution of Entropy
    - Distribution of Repetitiveness Score

    Each metric is plotted with transparent bars for better visualization of LHS and RHS
    distributions on the same axes. The x-axis limits for each histogram are manually set
    for specific ranges.

    Returns:
    - None

    Displays the histograms comparing LHS and RHS sequences for the calculated metrics.
    """
    
    max_homopolymer_lengths_lhs = lhs.apply(max_homopolymer_length)
    gc_contents_lhs = lhs.apply(calculate_gc_content)
    entropies_lhs = lhs.apply(calculate_entropy)
    repetitiveness_scores_lhs = lhs.apply(score_sequence_repetition)
    
    max_homopolymer_lengths_rhs = rhs.apply(max_homopolymer_length)
    gc_contents_rhs = rhs.apply(calculate_gc_content)
    entropies_rhs = rhs.apply(calculate_entropy)
    repetitiveness_scores_rhs = rhs.apply(score_sequence_repetition)
    
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    sns.histplot(max_homopolymer_lengths_lhs, bins=20, ax=axes[0, 0], color='darkorange', edgecolor='black', alpha=0.25, label='LHS')
    sns.histplot(max_homopolymer_lengths_rhs, bins=20, ax=axes[0, 0], color='cornflowerblue', edgecolor='black', alpha=0.25, label='RHS')
    axes[0, 0].set_xlabel('Max Homopolymer Length')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Max Homopolymer Length')
    axes[0, 0].legend()
    axes[0, 0].set_xlim(0, 6)  

    sns.histplot(gc_contents_lhs, bins=20, ax=axes[0, 1], color='darkorange', edgecolor='black', alpha=0.25, label='LHS')
    sns.histplot(gc_contents_rhs, bins=20, ax=axes[0, 1], color='cornflowerblue', edgecolor='black', alpha=0.25, label='RHS')
    axes[0, 1].set_xlabel('GC Content (%)')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of GC Content')
    axes[0, 1].legend()
    axes[0, 1].set_xlim(0, 100)  
    
    sns.histplot(entropies_lhs, bins=20, ax=axes[1, 0], color='darkorange', edgecolor='black', alpha=0.25, label='LHS')
    sns.histplot(entropies_rhs, bins=20, ax=axes[1, 0], color='cornflowerblue', edgecolor='black', alpha=0.25, label='RHS')
    axes[1, 0].set_xlabel('Entropy')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Entropy')
    axes[1, 0].legend()
    axes[1, 0].set_xlim(0.5, 1)  
    
    sns.histplot(repetitiveness_scores_lhs, bins=20, ax=axes[1, 1], color='darkorange', edgecolor='black', alpha=0.25, label='LHS')
    sns.histplot(repetitiveness_scores_rhs, bins=20, ax=axes[1, 1], color='cornflowerblue', edgecolor='black', alpha=0.25, label='RHS')
    axes[1, 1].set_xlabel('Repetitiveness Score')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Distribution of Repetitiveness Score')
    axes[1, 1].legend()
    axes[1, 1].set_xlim(0, 250)  

    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.show()