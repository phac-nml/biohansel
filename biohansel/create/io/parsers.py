from typing import Dict

import pandas as pd

from Bio import Seq, SeqIO

def parse_vcf(vcf_file: str) -> (pd.DataFrame, pd.DataFrame):
    """Reads in the generated vcf file, filters for 2-state SNVs and returns two dataframes: one that contains the
    REF/ALT sequence info and the other that contains just the binary SNV data for each genome
    Args: 
    vcf_file: path to the vcf file

    Returns:
    sequence_df: the DataFrame that contains just the REF/ALT sequence info
    binary_df: the DataFrame that contains the binary SNV data
    """
    with open(vcf_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#CHROM'):
                break
    cols = line.split()

    df = pd.read_table(
        vcf_file, comment='#', header=None,
        names=cols).rename(columns={'#CHROM': 'CHROM'})
    df = df[df['ALT'].str.len() <= 1]
    df = df.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], 1)
    df.columns = df.columns.str.replace("mysnps", "")
    df.index = df.POS
    sequence_df = df[['CHROM', 'REF', 'ALT']]
    binary_df = df.drop(['CHROM', 'POS', 'REF', 'ALT'], 1)

    return sequence_df, binary_df

def parse_sequence_file(reference_genome_path: str, reference_genome_type: str) -> Dict[str, Seq.Seq]:
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
