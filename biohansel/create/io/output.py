import os

import pandas as pd

from typing import List

from biohansel.create.schema_generator import get_sequence_string

def write_sequence_file(output_directory: str, df_list: List[pd.DataFrame],
                    schema_name: str, group: str) -> None:
    """Writes the sequences associated with each SNV into the output schema file
    Args:
        output_directory: directory where the schema would be located as indicated by the user
        df_list: list of DataFrames related for each particular chromosome included in the reference genome file 
        schema_name: the name of the output schema file
        group: the group membership of the current group of SNVs

    Returns:
         Creates schema file in the output directory specified by the user
    """

    with open(os.path.join(output_directory, f"{schema_name}.fasta"),
              "a+") as file:
        for curr_df in df_list:
            for index, row in curr_df.iterrows():
                ratio_value = row.iloc[3]
                position = index
                chromosome = row['CHROM']
                reference_snv = row['ref_sequences']
                alternate_snv = row['alt_sequences']
                sequence_string = get_sequence_string(ratio_value, chromosome, position, group, reference_snv,
                                                      alternate_snv)
                file.write(sequence_string)