import io
import os
import pandas as pd
from typing import List, Dict


def extract_test_columns(
        data_frame: pd.DataFrame, test_indices: list,
        reference_groups: Dict[str, str]) -> (pd.DataFrame, Dict[str, str]):
    """Removes any genomes that do not belong in the test group
    Args: 
    reference_groups: the groups file path that indicates the group that each genome belongs to
    data_frame: dataframe after being filtered for two-state SNVs

    Returns:
    new_data_frame: returns the genomes that are actually going to be used for the test
    test_dict: the dictionary that contains the group assignments for the 
    """
    test_dict = reference_groups
    data_frame.columns = data_frame.columns.str.strip()
    new_data_frame = data_frame

    for index in range(len(test_indices)):
        current_value = test_indices[index]
        new_data_frame = new_data_frame.drop(current_value, 1)
        test_dict.pop(current_value)

    return new_data_frame, test_dict
