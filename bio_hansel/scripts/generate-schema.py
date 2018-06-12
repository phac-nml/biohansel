import pandas as pd
import io
import os
from getsequence import getSequences


def main():
    with open('feht_resultsA.txt', 'r') as f:
        lines = [l for l in f if not l.startswith(('[#-', 'Group1', 'Group2', '--' , '-#]', '\n', 'Done'))]
        
       
         
       
    
    with open('feht_resultsA.txt', 'r') as file:
        all_lines=file.readlines()
    
    needed_string=all_lines[1].strip()
    character=needed_string[len(needed_string)-1:]


    print (character)
    data_frame=pd.read_table(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={'Name': int, 'GroupOne(+)': int, 'GroupOne(-)': int, 'GroupTwo(+)': int, 'GroupTwo(-)': int,
               'pValue': int, 'Ratio': int}
    )

    vcf=pd.read_table('filtered.vcf')

    new_data=data_frame.rename(columns={'Name':'POS'})
    result=new_data.merge(vcf[['POS','REF','ALT','CHROM']],how='left',on ='POS')
    getSequences(result, character)
      

           
        
   

if __name__ == '__main__':
    main()