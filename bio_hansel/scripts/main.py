import argparse
import logging
import os
import sys
import pandas as pd
import allel
# import numpy
# import scipy

SCRIPT_NAME='snv-subtyping-scheme'
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
    parser.add_argument('-i', '--input-vcf-file',
                        nargs=1,
                        metavar=('vcf-file'),
                        action='append',
                        help='vcf file of variants')
    parser.add_argument('-o', '--output-summary',
                        help='Subtyping summary output path (tab-delimited)')
    return parser
    
def main():
    parser = init_parser()
    # if len(sys.argv[1:]) == 0:
    #     parser.print_help()
    #     parser.exit()
    args = parser.parse_args()
    # init_console_logger(args.verbose)
    data_frame=allel.vcf_to_dataframe('core.vcf')
   
    data_frame=data_frame.query('ALT_1.len()>1')
    data_frame.to_csv('output.csv',sep='\t')
    print('completed action')


    # logging.debug(args)

if __name__ == '__main__':
    main()


    



