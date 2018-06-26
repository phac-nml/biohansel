from Bio import Entrez, SeqIO
import pandas as pd
import os


def getSequences(data_frame: pd.DataFrame, group:str, random_id:int, output_directory:str, reference_genome_path:str):
    """Collects the sequences from the from the reference genome by going through the dataframe and finding the appropriate SNV location
    Args: 
    output_directory:directory where the schema would be located as indicated by the user
    data_frame: filtered data frame with list of SNVs and their location
    group: specific group in which the SNV belongs to
    random_id: id that is assigned to the schema file
    output_directory: directory in which output files are stored
    reference_genome_path: file path to where the reference genome is located

    Output:
    Creates schema file in the output directory


    """
    max_sequence= data_frame.loc[data_frame['POS'].idxmax()]
    gb_file = reference_genome_path
    max_sequence_value=max_sequence['POS']
    with open(f"{output_directory}/schema-{random_id}.fasta", "a+") as file:
        for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
            for index, row in data_frame.iterrows():
                position=row['POS']
                reference_snv=row['REF']
                alternate_snv=row['ALT']

                seq_start1=max(0,position-16)-1
                seq_stop1=position-1

                seq_start2=position
                seq_stop2=min(max_sequence_value,position+16)

                record_1=gb_record.seq[seq_start1:seq_stop1]
                record_2=gb_record.seq[seq_start2:seq_stop2]
                #if the ratio is above 1, then it means that it is positive and takes the alternate snv form 
                if(row['Ratio']>0): 
                    file.write('>'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+alternate_snv+str(record_2)+'\n')
                    file.write('>negative'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+reference_snv+str(record_2)+'\n')

    
                    
                # if the ratio is below 1, then it means that it remains negative
                else:
                    file.write('>'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+reference_snv+str(record_2)+'\n')
                    file.write('>negative'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+alternate_snv+str(record_2)+'\n')

           

