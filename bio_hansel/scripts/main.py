import sys
import argparse
import logging
import wget
import os 
import pandas
from executeSnippy import executeSnippy
from filtervcf import filter_vcf, read_vcf, createSeparateVCF
from split_genomes import split
from getsequence import getSequences
from generateschema import generate_schema
from fishertest import conductFisherTest
from findcluster import findClusters


SCRIPT_NAME='schema_creation'
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)

def init_parser():
    parser = argparse.ArgumentParser(prog=SCRIPT_NAME,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                   
    parser.add_argument('-f', '--input_genomes',
                           required=True,
                        help='list of genomes to be analyzed by SRR accession (required)')
    parser.add_argument('-g', '--reference-genome-file',
                          required=True,
                        help='Path to reference genome name, can be .gb or .gbk format') 
    parser.add_argument('-o', '--output-folder-name',
                        help='Output file name for the genomes, default folder name is /usr/home/genomes')                   
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')
    return parser

def downloadFastqs(input_genomes: str, output_folder_name: str):
    home_folder = os.path.expanduser('~')
    """Downloads the FastQs from the given list of genomes from user input text file
       and outputs into user-specified output folder. If not, the FastQ files would go into a /genomes folder from the user's 
       home directory.

       Args:
       input_genomes: the path that leads to the genomes text file
       output_folder_name: the output folder name that the user has specified

       Returns:
       output_directory: the full output directory of where all of the output files will be stored
       genomes_file: The txt file with genomes that actually worked
    
    """
    output_directory=""
    try:
   
        with open (input_genomes, 'r') as file:
            print('Beginning file download with wget module')

            if output_folder_name is not None:
                
                output_directory=f"{home_folder}/{output_folder_name}"
            else:
                output_directory=f"{home_folder}/genomes"
                
            
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
            with open (f"{output_directory}/genomes_used.txt", "w") as genomes_file:
                for i in range(2):
                    while True:
                        try:  
                            for line in file:

                                print(f"output files will be stored in {output_directory}")
                                line=line.rstrip()
                                fastq1_string=f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{line[0:6]}/00{line[len(line)-1]}/{line}/{line}_1.fastq.gz"
                                fastq2_string=f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{line[0:6]}/00{line[len(line)-1]}/{line}/{line}_2.fastq.gz"
                                print(f"Reading in SRR name and attempting FastQ download for: {line}")
                                
                                file1=wget.download(fastq1_string, out=output_directory)
                                file2=wget.download(fastq2_string, out=output_directory)
                            if file1 is not None:
                                genomes_file.write(f"{line}\n")
                            return output_directory, genomes_file


            # handle error if files are not downloadable through wget
                        except IOError as err:
                            print("I/O error({0}): {1}".format(err.errno, err.strerror, input_genomes))
                        
                            print(f"There was a problem downloading the FastQ, {line} was not downloaded")
        
    #handle error for file input
    except IOError as e:
    # print(f"{e} was not found as a valid file name input")
        print("I/O error({0}): {1}".format(e.errno, e.strerror, input_genomes))
        print(e)
        sys.exit(1)
            

def main():
    home_folder = os.path.expanduser('~')
    parser = init_parser()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    output_folder_name= args.output_folder_name
    init_console_logger(3)
    input_genomes=args.input_genomes
    reference_genome_path=f"{home_folder}{args.reference_genome_file}"
   
    output_directory, input_genomes=downloadFastqs(input_genomes,output_folder_name)
    print(f"using genbank file from {reference_genome_path}")

    tree_file=executeSnippy(output_directory, reference_genome_path, input_genomes)
    groups_file=findClusters(tree_file, input_genomes, output_directory)
    data_frame=read_vcf(output_directory)
    test_indices=split(input_genomes)
    data_frame=filter_vcf(output_directory, data_frame)

    modified_data_frame, test_group=createSeparateVCF(data_frame, test_indices, output_directory, reference_groups, reference_genome_path)
    results_list=conductFisherTest(modified_data_frame, output_directory, test_group)
    generate_schema(output_directory, results_list, reference_genome_path)
    
if __name__ == '__main__':
    main()
