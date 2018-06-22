from Bio import Entrez, SeqIO
import pandas as pd
import os


def getSequences(data_frame, group, random_id, output_directory, reference_genome_path):
     """Generates the schema to be used as input for biohansel, goes through the results list and adds each group's individual SNVs into a schema file
    Args: 
    output_directory:directory where the schema would be located as indicated by the user
    results_list: the list of file paths that lead to the results of each group's fisher exact test ratio values
    reference_genome_path: the path to the genbank reference genome file

    """
    
    home_folder=os.path.expanduser('~')
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
                #     SeqIO.write(record, output_handle, "fasta")
                #if the reatio is above 1, then it means that it is positive and takes the alternate snv form 
                if(row['Ratio']>0): 
                    file.write('>'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+alternate_snv+str(record_2)+'\n')
                    file.write('>negative'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+reference_snv+str(record_2)+'\n')

                    # print('>'+str(position)+group+'\n')
                    # print(str(record_1)+alternate_snv+str(record_2)+'\n')
                    # print('>negative'+str(position)+group+'\n')
                    # print(str(record_1)+reference_snv+str(record_2)+'\n')
                    
                # if the ratio is below 1, then it means that it remains negative
                else:
                    file.write('>'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+reference_snv+str(record_2)+'\n')
                    file.write('>negative'+str(position)+'-'+group+'\n')
                    file.write(str(record_1)+alternate_snv+str(record_2)+'\n')

                    # print('>'+str(position)+'-'+group+'\n')
                    # print(str(record_1)+reference_snv+str(record_2)+'\n')
                    # print('>negative'+str(position)+'-'+group+'\n')
                    # print(str(record_1)+alternate_snv+str(record_2)+'\n')

