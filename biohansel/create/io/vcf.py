from typing import List, Tuple

import pandas as pd


def parse_vcf(vcf_file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Parse a VCF file into SNV info table and binary SNV matrix table

    SNV info table contains SNV state info (REF vs ALT) as well as position(POS) in the reference chromosome (CHROM).

    Binary SNV matrix contains SNV absence/presence in all samples.

    Notes:
        Both DataFrames have the same index.

    Args:
        vcf_file_path: VCF file path

    Returns:
        Tuple of
            - SNV state and location info
            - Binary SNV matrix; absence/presence of SNV in all samples
    """

    cols = _parse_vcf_headers(vcf_file_path)

    nt_columns = ['CHROM', 'POS', 'REF', 'ALT']
    assert len(set(nt_columns) - set(
        cols)) == 0, f'Not all expected nucleotide state information columns were parsed from {vcf_file_path}'
    unneeded_columns = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    assert len(set(unneeded_columns) - set(
        cols)) == 0, f'A standard VCF file should contain the following columns: {unneeded_columns}'

    df = pd.read_table(vcf_file_path, comment='#', header=None, names=cols)
    # keep only single nucleotide variations
    df = df[df['ALT'].str.len() == 1]
    # drop unneeded columns
    df = df.drop(unneeded_columns, 1)
    # split VCF table into table with SNV state info and table of the binary SNV matrix (presence/absence of SNV in all samples)
    df_nt = df[nt_columns]
    df_bin = df.drop(nt_columns, 1)
    # remove SNVs present in only one sample
    # row_sums = df_bin.sum(axis=1)
    # df_nt = df_nt.loc[row_sums > 1,]
    # df_bin = df_bin.loc[row_sums > 1,]
    return df_nt, df_bin


def _parse_vcf_headers(vcf_file_path: str) -> List[str]:
    header_line = None
    with open(vcf_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#CHROM'):
                header_line = line[1:]
                break
    if header_line is None:
        raise IOError(f'Could not find header line in VCF file "{vcf_file_path}"! Please ensure you have a correctly '
                      f'formatted VCF file.')
    cols = header_line.split()
    return cols
