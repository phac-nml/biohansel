import pandas as pd
import io
import os
from getsequence import getSequences
import uuid


def generate_schema(output_directory,results_list, reference_genome_path):
    random_id=uuid.uuid4()
    print(f"random if is {random_id}")

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
