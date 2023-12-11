import pandas as pd
from Bio.Seq import Seq
from datetime import datetime
import os

def generate_flex_sgrna_opool(a_position_protospacers, b_position_protospacers, pool_name, n_barcodes):
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
    sgA_RHS_CR = 'CTGAAAC'
    sgB_RHS_CR = 'AAAC'
    pCS1 = 'CGGTCCTAGCAA'

    ptruseq_r2 = "CAGACGTGTGCTCTTCCGATCT"
    lhs_cs = "ACGCGGTTAGCACGTANN"
    sgrna_a_cr = "CTTGCTATGCACTCTTGTGCTTAGCT"
    sgrna_b_cr = "GCTATGCTGTTTCCAGCTTAGCTCTT"

    a_seq = a_position_protospacers.apply(Seq)
    b_seq = b_position_protospacers.apply(Seq)

    pool_names_rhs = [f"RHS_{pool_name}"] * (len(a_seq) * 2)
    sequences_rhs = []

    for a, b in zip(a_seq, b_seq):
        sequences_rhs.extend([
            f"/5Phos/{sgA_RHS_CR}{a.reverse_complement()}{pCS1}",
            f"/5Phos/{sgB_RHS_CR}{b.reverse_complement()}{pCS1}"
        ])

    RHS_probe_order = pd.DataFrame({
        'Pool_name': pool_names_rhs,
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
            f"{ptruseq_r2}{lhs_cs}{barcode_sequences[i]}{sgrna_a_cr}",
            f"{ptruseq_r2}{lhs_cs}{barcode_sequences[i]}{sgrna_b_cr}"
        ])

    LHS_probe_order = pd.DataFrame({
        'Pool_name': pool_names_lhs,
        'Sequence': sequences_lhs
    })

    oligo_pool = pd.concat([RHS_probe_order, LHS_probe_order])

    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    directory = 'opools'
    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        oligo_pool.to_excel(writer, index=False)
