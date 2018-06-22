import pandas as pd
import io
import os
from getsequence import getSequences
import uuid


def generate_schema(output_directory: str,results_list: list, reference_genome_path: str):
     """Generates the schema to be used as input for biohansel, goes through the results list and adds each group's individual SNVs into a schema file
    Args: 
    output_directory:directory where the schema would be located as indicated by the user
    results_list: the list of file paths that lead to the results of each group's fisher exact test ratio values
    reference_genome_path: the path to the genbank reference genome file

    """
    random_id=uuid.uuid4()
    print(f"random id is {random_id}")

    for i in range(len(results_list)):
        with open(results_list[i], 'r') as file:
            lines = [line for line in file if not line.startswith(('[#-', 'Group1', 'Group2', '--' , '-#]', '\n', 'Done'))]
            

        with open(results_list[i], 'r') as file:
            all_lines=file.readlines()
        
        needed_string=all_lines[1].strip()
        character=needed_string[len(needed_string)-1:]


        print (character)
        data_frame=pd.read_table(
            io.StringIO(str.join(os.linesep, lines)),
            dtype={'Name': int, 'GroupOne(+)': int, 'GroupOne(-)': int, 'GroupTwo(+)': int, 'GroupTwo(-)': int,
                'pValue': int, 'Ratio': int}
        )

        vcf=pd.read_table(output_directory+'/filtered.vcf')

        new_data=data_frame.rename(columns={'Name':'POS'})
        result=new_data.merge(vcf[['POS','REF','ALT','CHROM']],how='left',on ='POS')
        getSequences(result, character, random_id, output_directory, reference_genome_path)
    print(f"job completed: output generated at schema-{random_id}.fasta")
