import os

from typing import Dict

import pandas as pd

from Bio import SeqIO


def get_sequences(results_dict: Dict[str, pd.DataFrame], sequence_length: int,
                ) -> Dict[str, pd.DataFrame]:
    """Collects the sequences from the reference genome by going through the list of DataFrames and adds two columns
    that contains the reference sequence and the alternate sequence surrounding each SNV

    Args:
        results_dict: specifies the list of genomes and their associated group
        sequence_length: the length of additional sequences to be added to the beginning and end of the SNV
        

    Returns:
        results_dict: updated dictionary with snv sequences for both the reference genome and alternate snv
    """
    for group, curr_df in results_dict.items():
        gb_record = [
            record for record in SeqIO.parse(reference_genome_path, "genbank")
        ]
        max_sequence_value = len(gb_record[0].seq)
        sequences = str(gb_record[0].seq)
        ref_seqs = curr_df.index.apply(
            get_sub_sequences,
            args=(sequences, sequence_length, max_sequence_value))
        alt_seqs = ref_seqs.str.slice(
            0, sequence_length) + curr_df.ALT + ref_seqs.str.slice(
                sequence_length + 1, sequence_length + sequence_length + 1)

        curr_df['ref_sequences'] = ref_seqs
        curr_df['alt_sequences'] = alt_seqs

    return results_dict


def write_sequences(output_directory: str,
                    updated_results_dict: Dict[str, pd.DataFrame],
                    schema_name: str) -> None:
    """Writes the sequences associated with each SNV into the output schema file
    Args:
        output_directory: directory where the schema would be located as indicated by the user
        updated_results_dict: updated dictionary with snv sequences for both the reference genome and alternate snv
        schema_name: the name of the output schema file

    Returns:
         Creates schema file in the output directory specified by the user
    """
    for group, curr_df in updated_results_dict.items():
        with open(os.path.join(output_directory,f"{schema_name}.fasta", "a+") as file:
            for i, row in curr_df.iterrows():
                attribute_value = row.iloc[2]
                position = i
                reference_snv = row['ref_sequences']
                alternate_snv = row['alt_sequences']
                # if the ratio is above 1, then it means that it is positive and takes the alternate snv form
                if attribute_value > 0:
                    file.write(f'''
                            >{position}-{group}
                            {alternate_snv}
                            >negative{position}-{group}
                            {reference_snv}
                            ''')
                # if the ratio is below 1, then it means that it remains negative
                else:
                    file.write(f'''
                            >{position}-{group}
                            {reference_snv}
                            >negative{position}-{group}
                            {alternate_snv}
                            ''')


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
