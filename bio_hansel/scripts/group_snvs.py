from typing import Dict

import pandas as pd


def group_snvs(binary_df: pd.DataFrame, sequence_df: pd.DataFrame
                test_groups: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """Takes in a DataFrame containing SNV VCF data and extracts the SNVs that are specific to a group, and only to that
    group

    Args:
        binary_df: the DataFrame that contains the binary SNV data
        sequence_df: the DataFrame that contains just the REF/ALT sequence info
        groups_dict: the dictionary that contains the group information for each genome
    Returns:
        results_list: A dictionary containing the group allocation and DataFrame of SNVs that are associated with that
                      group
    """
   
    unique_groups = list(set(test_groups.values()))
    results_list = {}
    other_list = []
    current_list = []
    for x in unique_groups:
        for key, value in test_groups.items():
            if x == value:
                current_list.append(key)
            else:
                other_list.append(key)
        dfsnv_curr = binary_df[current_list]
        dfsnv_other = binary_df[other_list]
        row_sums_curr = dfsnv_curr.sum(axis=1)
        row_sums_other = dfsnv_other.sum(axis=1)
        new_data_frame = (dfsnv_curr.loc[(
                                                 (row_sums_curr == 0) & (row_sums_other == len(other_list))
                                         ) | ((row_sums_curr == len(current_list)) & (row_sums_other == 0)), :])
        final_table = pd.concat([sequence_df, new_data_frame], axis=1)
        final_table = final_table[final_table.columns[:3]]
        results_list[x] = final_table.dropna()
        current_list = []
        other_list = []

    return results_list
