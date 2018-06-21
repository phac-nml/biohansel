import os
from subprocess import call
import shutil


def executeSnippy(output_directory: str, reference_genome: str, inputgenomes: str):
    """Runs the external Snippy command-line tool that will take the downloaded genomes and creat

    """
    snippy_command="snippy1"
    snippycore_command="snippy-core"
    try:
        if(shutil.which(snippy_command)) is not None:
        
            with open(inputgenomes,'r') as file:
                for line in file:
                    line=line.rstrip()
                    call(["snippy", "--outdir",f"{output_directory}/mysnps{line}", "--ref",reference_genome, "--R1",f"{output_directory}/{line}_1.fastq.gz", "--R2", f"{output_directory}/{line}_2.fastq.gz"]) 
                
                argument_list=["snippy-core", "--prefix", f"{output_directory}/core"]

                for line in file:
                    line=line.rstrip()
                    argument_list.append(f"{output_directory}/mysnps{line}")
                    
                call(argument_list) 
        else:
            raise Exception("snippy not installed on this machine, cannot proceed any further")
            
    except Exception as e:
        print(e.args)

    return (f"{output_directory}/referencegroups.txt")

        ##still need to use FastTree and then create a testreference.txt file
    
