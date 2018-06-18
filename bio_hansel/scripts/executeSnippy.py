import os
from subprocess import call


def executeSnippy(output_directory, reference_genome, inputgenomes):

    
    with open(inputgenomes,'r') as file:
    #     for line in file:
    #         line=line.rstrip()
    #         call(["snippy", "--outdir",output_directory+"/mysnps"+line, "--ref",reference_genome, "--R1",output_directory+"/"+line+"_1.fastq.gz", "--R2", output_directory+"/"+line+"_2.fastq.gz"]) 
        

        argument_list=["snippy-core", "--prefix", output_directory+"/core"]

        for line in file:
            line=line.rstrip()
            argument_list.append(output_directory+"/mysnps"+line)
            


        call(argument_list) 

        return (f"{output_directory}/referencegroups.txt")

        ##still need to use FastTree and then create a testreference.txt file
    

if __name__ == '__main__':
    executeSnippy()
