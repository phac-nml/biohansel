import math
import os


from typing import Dict, List

import pandas as pd

from Bio import SeqIO, Seq


def get_sequences(
        curr_df: pd.DataFrame,
        tile_length: int,
        record_dict: Dict[str, Seq.Seq],
) -> List[pd.DataFrame]:
    """Collects the sequences from the reference genome by going through the list of DataFrames and adds two columns
    that contains the reference sequence and the alternate sequence surrounding each SNV

    Args:
        curr_df: a DataFrame that contains each group's individual SNV data
        sequence_length: the length of additional sequences to be added to the beginning and end of the SNV
        record_dict: dictionary of records obtained from the reference sequence file

    Returns:
        df_list: list of DataFrames that contain snv sequences related to each chromosome for that particular cluster group
    """
    pad_length=math.ceil((tile_length-1)/2)
    df_list = []
    curr_df.CHROM = curr_df.CHROM.str.strip()
    unique_chromosomes = curr_df.CHROM.unique()
    for chromosome in unique_chromosomes:
        record_seq = record_dict[chromosome]
        max_sequence_value = len(record_seq)
        sequences = str(record_seq)
        curr_chrom_df = curr_df[curr_df['CHROM'].str.match(
            chromosome.strip())]

        positions = pd.Series(curr_chrom_df.index, index=curr_chrom_df.index)
        ref_seqs = positions.apply(
            get_subsequences,
            args=(sequences, pad_length, max_sequence_value))
        alt_seqs = ref_seqs.str.slice(
            0, pad_length) + curr_chrom_df.ALT + ref_seqs.str.slice(
            pad_length + 1, pad_length + pad_length + 1)

        curr_chrom_df['ref_sequences'] = ref_seqs
        curr_chrom_df['alt_sequences'] = alt_seqs
        df_list.append(curr_chrom_df)

    return df_list

def get_subsequences(position: int, seq: str, pad_length: int,
                     max_sequence_value: int) -> str:
    """Get the sequences that are before are after the specified SNV
    Args:
        position: the position of the SNV based off of the reference genome
        seq: the set of sequences from the reference genome
        sequence_length: the amount of upstream and downstream sequences added to the SNV and included in the schema
                        output
        max_sequence_value: the length of the reference genome

    Returns:
        sequence_sequence: the sequence that is found at the specified genome positions from the reference genome

    """

    specific_sequence = seq[max(0, position - (pad_length + 1)):min(
        max_sequence_value, position + pad_length)]

    return specific_sequence


def get_sequence_string(ratio_value: int, chromosome: str, position: int, group: str, reference_snv: str,
                        alternate_snv: str) -> str:
    """Creates the string output for the schema file
    Args:
        ratio_value: the ratio_value for that particular snv, with 1 being that all genomes contain the snv and 0 being
                    none of the genomes contain the sequence
        chromosome: the current chromosome of the reference genome sequence
        position: the position of the SNV
        group: the group membership of the SNV
        reference_snv: the genome sequence of the reference genome
        alternate_snv: the genome sequence of the alternate form of the SNV

    Returns:
        sequence_string: the sequence string to be outputted into the schema file

    """


    if (len(str(group))>1):
        sequence_string=f""">{position}-{group}
{alternate_snv if ratio_value > 0 else reference_snv}
>negative{position}-{group}
{reference_snv if ratio_value > 0 else alternate_snv}
"""
    else:
        sequence_string=f""">{position}-{group}
{alternate_snv if ratio_value > 0 else reference_snv}
"""
    return sequence_string


    


def group_snvs(
        binary_df: pd.DataFrame,
        sequence_df: pd.DataFrame,
        groups_dict: Dict[str, str],
) -> Dict[str, pd.DataFrame]:
    """Takes in a DataFrame containing SNV VCF data and extracts the SNVs that are specific to a group, and only to that
    group

    Args:
        binary_df: the DataFrame that contains the binary SNV data
        sequence_df: the DataFrame that contains just the REF/ALT sequence info
        groups_dict: the dictionary that contains the group information for each genome
    Returns:
        results_list: A dictionary containing the group allocation and DataFrame of SNVs that are associated with that
                      group
    """

    unique_groups = list(set(groups_dict.values()))
    results_list = {}
    other_list = []
    current_list = []

    for group in unique_groups:
        for genome, curr_group in groups_dict.items():
            if group == curr_group:
                current_list.append(genome)
            else:
                other_list.append(genome)
        dfsnv_curr = binary_df[current_list]
        dfsnv_other = binary_df[other_list]
        row_sums_curr = dfsnv_curr.sum(axis=1)
        row_sums_other = dfsnv_other.sum(axis=1)
        distinct = (row_sums_curr == 0) & (row_sums_other == len(other_list))
        all_negative = (row_sums_curr == len(current_list)) & (row_sums_other == 0)
        group_snv_df = dfsnv_curr.loc[distinct | all_negative, :]
        final_table = pd.concat([sequence_df, group_snv_df], axis=1)
        final_table = final_table[final_table.columns[:4]]
        results_list[group] = final_table.dropna()
        current_list = []
        other_list = []

    return results_list
