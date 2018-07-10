import io
import os
import pandas as pd
from typing import List, Dict
import logging


def extract_test_columns(
        df: pd.DataFrame, test_indices: list,
        reference_groups: Dict[str, str]) -> (pd.DataFrame, Dict[str, str]):
    """Removes any genomes that do not belong in the test group according to the test_indices
    Args: 
    reference_groups: the groups file path that indicates the group that each genome belongs to
    df: dataframe after being filtered for two-state SNVs

    Returns:
    new_df: returns the genomes that are actually going to be used for the test
    test_dict: the dictionary that contains the group assignments for the 
    """
    test_dict = reference_groups
    df.columns = df.columns.str.strip()
    new_df = df
    logging.debug('Reference groups: %s', reference_groups)

    for index in range(len(test_indices)):
        current_value = test_indices[index]
        new_data_frame = new_df.drop(current_value, 1)
        test_dict.pop(current_value)

    return new_data_frame, test_dict
