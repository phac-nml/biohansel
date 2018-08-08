import os
import textwrap

from typing import Dict, List

import pandas as pd

from Bio import SeqIO, Seq


def get_sequences(
        curr_df: pd.DataFrame,
        sequence_length: int,
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
        group: the group membership of the current group of SNVs

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
                sequence_string = get_sequence_string(ratio_value, chromosome, position, group, reference_snv,
                                                      alternate_snv)
                file.write(sequence_string)


def get_subsequences(position: int, seq: str, sequence_length: int,
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


def read_sequence_file(reference_genome_path: str, reference_genome_type: str) -> Dict[str, Seq.Seq]:
    """Reads in the sequence file and indexes each of the individual sequences into a dictionary 
    to allow for faster querying
    Args:   
        reference_genome_path: the path to the reference genome
        reference_genome_type: reference genome file type

    Returns:
        record_dict: returns a dictionary of all the record sequences indexed by record name
    """
    record_dict = {}
    for record in SeqIO.parse(reference_genome_path, reference_genome_type):
        record_dict[record.name] = record.seq

    return record_dict


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

    return (textwrap.dedent(f"""\
    >({chromosome}){position}-{group}
    {alternate_snv if ratio_value > 0 else reference_snv}
    >negative({chromosome}){position}-{group}
    {reference_snv if ratio_value > 0 else alternate_snv}\n"""))
