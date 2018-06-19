from Bio import Entrez, SeqIO
import pandas as pd
import os


def getSequences(data_frame, group, random_id, output_directory, reference_genome_path):
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
                chromosome=row['CHROM']
            
                # Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
                # handle_1 = Entrez.efetch(db="nucleotide", 
                #                     id=chromosome, 
                #                     rettype="fasta", 
                #                     strand=1, 
                #                     seq_start=max(0,position-17), 
                #                     seq_stop=position-1)
                

                # handle_2 = Entrez.efetch(db="nucleotide", 
                #                     id=chromosome, 
                #                     rettype="fasta", 
                #                     strand=1, 
                #                     seq_start=position+1,
                #                     seq_stop=min(max_sequence_value,position+17))
                # record_1 = SeqIO.read(handle_1, "fasta")
                # record_2=SeqIO.read(handle_2,"fasta")

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
   
                
        
    
    


if __name__ == '__main__':
    getSequences()
