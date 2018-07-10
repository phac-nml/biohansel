import sys
import argparse
import logging
import wget
import os
from read_vcf import read_vcf
from extract_test_columns import extract_test_columns
from split_genomes import split_genomes
from get_sequence import get_sequences
from fisher_test import fisher_test
from find_cluster import find_clusters

SCRIPT_NAME = 'schema_creation'
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)


def init_parser():
    parser = argparse.ArgumentParser(
        prog=SCRIPT_NAME, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        '-i',
        '--input-vcf',
        required=True,
        help=
        'vcf file containing the list of snps that are related to the genomes')

    parser.add_argument(
        '-g',
        '--reference-genome-file',
        required=True,
        help='Path to reference genome name, can be .gb or .gbk format')
    parser.add_argument(
        '-o',
        '--output-folder-name',
        required=True,
        help=
        'Output file name for the genomes, default folder name is /usr/home/genomes'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help=
        'Logging verbosity level (-v == show warnings; -vvv == show debug info)'
    )
    return parser


def main():
    home_folder = os.path.expanduser('~')
    parser = init_parser()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    output_folder_name = args.output_folder_name
    init_console_logger(3)
    vcf_file = args.input_vcf
    reference_genome_path = f"{home_folder}{args.reference_genome_file}"
    print(f"using genbank file from {reference_genome_path}")

    output_directory = f"{home_folder}/{output_folder_name}"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    data_frame = read_vcf(vcf_file)
    temporary = data_frame.drop(['POS', 'REF', 'ALT'], 1)
    genomes_only = temporary.columns
    groups_dict = find_clusters(data_frame)

    test_indices = split_genomes(genomes_only)

    # this function basically takes in the data_frames and returns a dataframe with only the test genomes
    modified_data_frame, test_group = extract_test_columns(
        data_frame, test_indices, groups_dict)

    ## need to extract the all-or-nothing snvs and also to parse out the columns that are in the test groups
    results_dict = fisherTest(modified_data_frame, test_group)
    get_sequences(output_directory, reference_genome_path, results_dict)


if __name__ == '__main__':
    main()
