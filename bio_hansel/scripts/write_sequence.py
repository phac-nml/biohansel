import os

from typing import Dict, List

import pandas as pd

from Bio import SeqIO


def get_sequences(
        curr_df: pd.DataFrame,
        sequence_length: int,
        record_dict: Dict[str, str],
) -> List[pd.DataFrame]:
    """Collects the sequences from the reference genome by going through the list of DataFrames and adds two columns
    that contains the reference sequence and the alternate sequence surrounding each SNV

    Args:
        curr_df: a DataFrame that contains each group's individual SNV data
        sequence_length: the length of additional sequences to be added to the beginning and end of the SNV
        record_dict: dictionary of records obainted from the reference sequence file
        

    Returns:
        resis ults_dict: updated dictionary with snv sequences for both the reference genome and alternate snv
    """
    df_list = []
    curr_df.CHROM=curr_df.CHROM.str.strip()
    unique_chromosomes=curr_df.CHROM.unique()
    for chromosome in unique_chromosomes:
        record_seq=record_dict[chromosome]
        max_sequence_value = len(record_seq)
        sequences = str(record_seq)
        curr_chrom_df = curr_df[curr_df['CHROM'].str.match(
            chromosome.strip())]
        curr_chrom_df['POS'] = curr_chrom_df.index
        ref_seqs = curr_chrom_df.POS.apply(
            get_sub_sequences,
            args=(sequences, sequence_length, max_sequence_value))
        alt_seqs = ref_seqs.str.slice(
            0, sequence_length) + curr_chrom_df.ALT + ref_seqs.str.slice(
                sequence_length + 1, sequence_length + sequence_length + 1)

        curr_chrom_df['ref_sequences'] = ref_seqs
        curr_chrom_df['alt_sequences'] = alt_seqs
        df_list.append(curr_chrom_df)
    return df_list


def write_sequences(output_directory: str, df_list: List[pd.DataFrame],
                    schema_name: str, group: str) -> None:
    """Writes the sequences associated with each SNV into the output schema file
    Args:
        output_directory: directory where the schema would be located as indicated by the user
        df_list: list of DataFrames related for each particular chromosome included in the reference genome file 
        schema_name: the name of the output schema file

    Returns:
         Creates schema file in the output directory specified by the user
    """

    with open(os.path.join(output_directory, f"{schema_name}.fasta"),
              "a+") as file:
        for curr_df in df_list:
            for index, row in curr_df.iterrows():
                ratio_value = row.iloc[3]
                position = index
                chromosome = row['CHROM']
                reference_snv = row['ref_sequences']
                alternate_snv = row['alt_sequences']
                # if the ratio is above 1, then it means that it is positive and takes the alternate snv form
                if ratio_value > 0:
                    file.write(
                        f">({chromosome}){position}-{group}\n"
                        f"{alternate_snv}\n"
                        f">negative({chromosome}){position}-{group}\n"
                        f"{reference_snv}\n")
                # if the ratio is below 1, then it means that it remains negative
                else:
                    file.write(
                        f">({chromosome}){position}-{group}\n"
                        f"{reference_snv}\n"
                        f">negative({chromosome}){position}-{group}\n"
                        f"{alternate_snv}\n")


def get_sub_sequences(position: int, seq: str, sequence_length: int,
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
    specific_sequence = seq[max(0, position - (sequence_length + 1)):min(
        max_sequence_value, position + sequence_length)]
    return specific_sequence

def read_sequence_file(reference_genome_path: str)-> Dict[str, str]:
    """Reads in the sequence file and indexes each of the individual sequences into a dictionary 
    to allow for faster querying
    Args:   
        reference_genome_path: the path to the reference genome

    Returns:
        record_dict: returns a dictionary of all the record sequences indexed by record name
        

    """
    record_dict={}
    for record in SeqIO.parse(reference_genome_path, "genbank"):
        
        record_dict[record.name]=record.seq
    
    return record_dict



