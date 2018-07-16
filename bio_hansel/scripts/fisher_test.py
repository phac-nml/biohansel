from typing import Dict

import pandas as pd


def fisher_test(modified_df: pd.DataFrame,
                test_groups: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """Takes in a DataFrame containing SNV VCF data and extracts the SNVs that are specific to a group, and only to that
    group

    Args:
        modified_df: a DataFrame containing SNV VCF data
        test_groups: the group of SNVs that are to be included for testing

    Returns:
        results_list: A dictionary containing the group allocation and DataFrame of SNVs that are associated with that
                      group
    """
    attributes = modified_df[['POS', 'REF', 'ALT']]
    snvs_only = modified_df.drop(['POS', 'REF', 'ALT'], axis=1)
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
        dfsnv_curr = snvs_only[current_list]
        dfsnv_other = snvs_only[other_list]
        row_sums_curr = dfsnv_curr.sum(axis=1)
        row_sums_other = dfsnv_other.sum(axis=1)
        new_data_frame = (dfsnv_curr.loc[(
                                                 (row_sums_curr == 0) & (row_sums_other == len(other_list))
                                         ) | ((row_sums_curr == len(current_list)) & (row_sums_other == 0)), :])
        final_table = pd.concat([attributes, new_data_frame], axis=1)
        final_table = final_table[final_table.columns[:4]]
        results_list[x] = final_table.dropna()
        current_list = []
        other_list = []

    return results_list
