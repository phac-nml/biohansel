import pandas as pd


def read_vcf(vcf_file: str) -> (pd.DataFrame):
    """Reads in the generated vcf file, filters for 2-state SNVs and takes out unnecessary columns within the vcf file
    Args: 
    vcf_file: path to the vcf file

    Returns:
    df: returns the dataframe version of the generated vcf file from snippy
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
    df = df.drop(['CHROM', 'ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], 1)

    return df
