import logging
import math
import os


from typing import Dict, List, Set

import pandas as pd

from Bio import SeqIO, Seq


def get_sequences(
        curr_df: pd.DataFrame,
        tile_length: int,
        record_dict: Dict[str, Seq.Seq],
) -> List[pd.DataFrame]:
    """Collects the sequences from the reference genome by going through the list of DataFrames and adds two columns
    that contains the reference sequence and the alternate sequence surrounding each snv

    Args:
        curr_df: a DataFrame that contains each group's individual snv data
        sequence_length: the length of additional sequences to be added to the beginning and end of the snv
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
    """Get the sequences that are before are after the specified snv
    Args:
        position: the position of the snv based off of the reference genome
        seq: the set of sequences from the reference genome
        sequence_length: the amount of upstream and downstream sequences added to the snv and included in the schema
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
        position: the position of the snv
        group: the group membership of the snv
        reference_snv: the genome sequence of the reference genome
        alternate_snv: the genome sequence of the alternate form of the snv

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

def find_snvs_discriminating_2_groups(binary_df: pd.DataFrame, ingroup_genomes: List[str], outgroup_genomes: List[str]) -> List[int]: 
    """Finds snvs that discriminates between the in-group genomes and outgroup genomes
    Args:
        binary_df:the DataFrame that contains the binary snv data
        ingroup_genomes: the list of in-group genomes
        outgroup_genomes: the list of out-group genomes

    Returns:
        distinct_list_of_snvs list of snv positions that separate the 2 groups of genomes
    
    """
    
    snv_df_ingroup = binary_df[ingroup_genomes]
    snv_df_outgroup = binary_df[outgroup_genomes]
    row_sums_ingroup = snv_df_ingroup.sum(axis=1)
    row_sums_outgroup = snv_df_outgroup.sum(axis=1)
    distinct_snvs = (row_sums_ingroup == 0) & (row_sums_outgroup == len(outgroup_genomes))
    all_negative_snvs = (row_sums_ingroup == len(ingroup_genomes)) & (row_sums_outgroup == 0)
    all_negative_set=set(all_negative_snvs.where(all_negative_snvs==True).dropna().index)
    all_distinct_set=set(distinct_snvs.where(distinct_snvs==True).dropna().index)
    distinct_list_of_snvs=list(all_negative_set.union(all_distinct_set))


    return distinct_list_of_snvs
def group_snvs(
        binary_df: pd.DataFrame,
        sequence_df: pd.DataFrame,
        groups_dict: Dict[str, str],
) -> Dict[str, pd.DataFrame]:
    """Takes in a DataFrame containing snv VCF data and extracts the snvs that are specific to a group, and only to that
    group

    Args:
        binary_df: the DataFrame that contains the binary snv data
        sequence_df: the DataFrame that contains just the REF/ALT sequence info
        groups_dict: the dictionary that contains the group information for each genome
    Returns:
        snvs_dict: A dictionary containing the group allocation and DataFrame of snvs that are associated with that
                      group
    """

    unique_groups = list(set(groups_dict.values()))
    snvs_dict = {}
    outgroup_genomes = []
    ingroup_genomes = []

    for group in unique_groups:
        for genome, curr_group in groups_dict.items():
            if group == curr_group:
                ingroup_genomes.append(genome)
            else:
                outgroup_genomes.append(genome)
        list_of_discriminating_snvs=find_snvs_discriminating_2_groups(binary_df, ingroup_genomes, outgroup_genomes)
        snv_df_ingroup=binary_df[ingroup_genomes]       
        group_snv_df = snv_df_ingroup.loc[pd.Series(list_of_discriminating_snvs), :]
       
        distinct_snv_df = pd.concat([sequence_df, group_snv_df], axis=1)
        distinct_snv_df = distinct_snv_df[distinct_snv_df.columns[:4]]
        distinct_snv_df.columns=['CHROM', 'REF', 'ALT', 'ratio_value']
        snvs_dict[group] = distinct_snv_df.dropna()
        ingroup_genomes = []
        outgroup_genomes = []
   
    return snvs_dict

