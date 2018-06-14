import sys
import argparse
import logging
import wget
import os 

SCRIPT_NAME='SRR_name_conversion'
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
                                    #  description=program_desc)
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

def main():
    parser = init_parser()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    init_console_logger(3)
    filename=args.input_genomes
    print(filename)
    file= open(filename,'r')
    

    print('Beginning file download with wget module')

    ##import files into home directory
    home_directory=os.path.expanduser('~')
    if args.output_folder_name is not None:
        
        output_directory=home_directory+'/'+args.output_folder_name
    else:
        output_directory=home_directory+'/'+'genomes'

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)


    for line in file:
        line=line.rstrip()
        newString1='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+line[0:6]+'/00'+line[len(line)-1]+'/'+line+'/'+line+'_1.fastq.gz '
        newString2='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+line[0:6]+'/00'+line[len(line)-1]+'/'+line+'/'+line+'_2.fastq.gz'
        wget.download(newString1, out=output_directory)
        wget.download(newString2, out=output_directory)
        


    file.close()
    

if __name__ == '__main__':
    main()
