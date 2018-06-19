import io
import os
import pandas as pd
from split_genomes import split
from fishertest import conductFisherTest
from generateschema import generate_schema



def read_vcf(output_directory):
    vcf_file=output_directory+"/core.vcf"
    with open(vcf_file, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]
    return pd.read_table(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str}
    ).rename(columns={'#CHROM': 'CHROM'})

def filter_vcf(output_directory, input_genomes, reference_groups, reference_genome_path):
    data_frame=read_vcf(output_directory)
    # print(data_frame)
    ##take out any columns that contains more than one state
    data_frame=data_frame[data_frame['ALT'].str.len()<=1]
    # print(data_frame)

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
    train_indices, test_indices=split(input_genomes)
    print(train_indices)

    createSeparateVCF(data_frame, test_indices, output_directory, reference_groups, reference_genome_path)


def createSeparateVCF(data_frame, test_indices, output_directory, reference_groups, reference_genome_path):
    new_data_frame=data_frame
    
    test_group=pd.read_table(reference_groups, sep='\t')
    test_group.columns=test_group.columns.str.strip()
    for index in range(len(test_indices)):
        
        current_value='mysnps'+test_indices[index]
        new_data_frame=new_data_frame.drop(current_value, 1)
        test_group=test_group[test_group.genomes!=current_value]

    # print(list(test_group.columns.values))    
    print(new_data_frame)  
    print(test_group)  
    
   
    modified_df=new_data_frame.drop(['CHROM','ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], 1)
    results_list=conductFisherTest(modified_df, output_directory, test_group)
    generate_schema(output_directory, results_list, reference_genome_path)

    




if __name__ == '__main__':
    main()