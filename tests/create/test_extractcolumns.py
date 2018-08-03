# import pandas as pd
# from bio_hansel.scripts.extract_test_columns import extract_test_columns
# from pandas.util.testing import assert_frame_equal
# test_data = pd.read_csv('tests/data/expected.csv')
# test_indices = ['mysnpsSRR6683541', 'mysnpsSRR6683736']
# test_dict = {
#     'Reference': 'A',
#     'mysnpsSRR6683541': 'B',
#     'mysnpsSRR6683736': 'D',
#     'mysnpsSRR6683914': 'A',
#     'mysnpsSRR6683916': 'B'
# }
# expected_group = pd.read_csv('tests/data/expected_aftercut.csv')
# expected_dict = {
#     'Reference': 'A',
#     'mysnpsSRR6683914': 'A',
#     'mysnpsSRR6683916': 'B'
# }


# def test_extract_columns():
#     resulting_data_frame, result_test_groups = extract_test_columns(
#         test_data, test_indices, test_dict)
#     assert_frame_equal(resulting_data_frame, expected_group)

#     assert (len(expected_dict) == len(result_test_groups))
#     assert (expected_dict == result_test_groups)
