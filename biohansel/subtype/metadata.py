import logging
from typing import Optional

import numpy as np
import os
import pandas as pd


def read_metadata_table(path: str) -> Optional[pd.DataFrame]:
    """Read a subtype metadata table into a Pandas DataFrame.

    This table must have a column labeled `subtype` and have one of the following file extensions:

    - `.tab` - tab-delimited file
    - `.tsv` - tab-delimited file
    - `.csv` - comma-separated values file

    The top row must be the header row.

    Args:
        path: File path of table.

    Returns:
        DataFrame of table file if `path` is one of the acceptable file formats, otherwise, return `None`
    """
    FILE_EXT_TO_PD_READ_FUNC = {
        '.tab': pd.read_table,
        '.tsv': pd.read_table,
        '.csv': pd.read_csv
    }
    _, file_ext = os.path.splitext(os.path.basename(path))
    file_ext = file_ext.lower()
    if file_ext not in FILE_EXT_TO_PD_READ_FUNC:
        logging.error('File extension of metadata file "{}" not one of the expected "{}"'.format(
            path,
            list(FILE_EXT_TO_PD_READ_FUNC.keys())
        ))
        return None
    dfmd = FILE_EXT_TO_PD_READ_FUNC[file_ext](path)  # type: pd.DataFrame
    assert np.any(dfmd.columns == 'subtype'), 'Column with name "subtype" expected in metadata file "{}"'.format(path)
    dfmd.subtype = dfmd.subtype.astype(str)
    logging.info('Read scheme metadata file "{}" into DataFrame with shape {}'.format(path, dfmd.shape))
    return dfmd


def merge_metadata_with_summary_results(df_results: pd.DataFrame, df_metadata: pd.DataFrame) -> pd.DataFrame:
    """Merge subtype metadata table into subtype results table.

    Args:
        df_results: Subtyping results table.
        df_metadata: Subtype metadata table.

    Returns:
        Subtyping results with subtype metadata merged in if metadata is present for subtype results.
    """
    df_results.subtype = df_results.subtype.fillna('')
    df_results.subtype = df_results.subtype.astype(str)
    return pd.merge(df_results, df_metadata, how='left', on='subtype')
