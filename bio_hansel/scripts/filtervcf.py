import io
import os
import pandas as pd
from split_genomes import split
from fishertest import conductFisherTest
from generateschema import generate_schema



def read_vcf(output_directory: str):
    """Reads in the generated vcf file
    Args: 
    Output_directory: where the output files are stored

    Returns:
    data_frame: returns the dataframe version of the generated vcf file from snippy
    """
    vcf_file=output_directory+"/core.vcf"
    with open(vcf_file, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]
    data_frame=pd.read_table(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str}
    ).rename(columns={'#CHROM': 'CHROM'})
    return data_frame

def filter_vcf(output_directory: str, data_frame:pd.DataFrame):
    """Removes any SNVs that have more than than two states, i.e. A, T, G
    Args: 
    Output_directory: where the output files are stored
    data_frame: A copy of the data_frame after it has been read in from core.vcf

    Returns:
    data_frame: returns the dataframe after filtering only for two-state SNVs


    """
    ##take out any columns that contains more than one state
    data_frame=data_frame[data_frame['ALT'].str.len()<=1]
    

    header = """##fileformat=VCFv4.1
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##source=myImputationProgramV3.1
    ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele.">
    #CHROM POS ID REF ALT QUAL FILTER INFO
    """

    output_VCF = output_directory+"/filtered.vcf"
    with open(output_VCF, 'w') as vcf:
      vcf.write(header)
    
    data_frame.to_csv(output_VCF, sep="\t", mode='w', index=False)
    return data_frame
    


def createSeparateVCF(data_frame: pd.DataFrame, test_indices: list, reference_groups: str):
    """Removes any SNVs that have more than than two states, i.e. A, T, G
    Args: 
    reference_groups: the groups file path that indicates the group that each genome belongs to
    data_frame: dataframe after filtering only for two-state SNVs

    Returns:
    data_frame: returns the genomes that are actually going to be used for the test


    """

    new_data_frame=data_frame

    test_group=pd.read_table(reference_groups, sep='\t')
    test_group.columns=test_group.columns.str.strip()
    for index in range(len(test_indices)):

        current_value='mysnps'+test_indices[index]
        new_data_frame=new_data_frame.drop(current_value, 1)
        test_group=test_group[test_group.genomes!=current_value]

    modified_df=new_data_frame.drop(['CHROM','ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], 1)
    return modified_df, test_group

