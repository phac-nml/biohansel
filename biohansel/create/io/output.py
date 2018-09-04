import os

import logging
import pandas as pd

from typing import List, Iterable

from biohansel.create.schema_generator import get_sequence_string

def write_sequence_file(output_directory: str,
                    schema_name: str, fasta_sequences: Iterable[str]) -> None:
    """Writes the sequences associated with each SNV into the output schema file
    Args:
        output_directory: directory where the schema would be located as indicated by the user
        schema_name: the name of the output schema file
        fasta_sequences: list of sequences to be added to the tiles file

    Returns:
         Creates schema file in the output directory specified by the user
    """

    with open(os.path.join(output_directory, f"{schema_name}.fasta"),
              "a+") as file_out:
        for sequence_string in fasta_sequences:
                file_out.write(sequence_string)

def generate_tile_fasta(df_list: List[pd.DataFrame], group: str) -> Iterable[str]:
    """Generates the list of fasta sequences associated with that DataFrame
    Args:
        df_list: list of DataFrames related for each particular chromosome included in the reference genome file 
        group: the group membership of the current group of SNVs
    Returns:
        list of fasta sequences
    """
    for curr_df in df_list:
        for index, row in curr_df.iterrows():
            ratio_value = row['ratio_value']
            position = index
            chromosome = row['CHROM']
            reference_snv = row['ref_sequences']
            alternate_snv = row['alt_sequences']
            yield get_sequence_string(ratio_value, chromosome, position, group, reference_snv, alternate_snv)
