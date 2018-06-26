import pandas as pd
from subprocess import call
import shutil


def conductFisherTest(modified_data_frame: pd.DataFrame, output_directory: str, test_groups:list):
    """Conducts Fisher's exact test on the dataframe, uses the external library: feht (https://github.com/chadlaing/feht)
    Feht must be previously downloaded on the user's machine
    Args: 
    modified_data_frame: the modified dataframe that has been filtered for SNVs 

    Returns:
    data_frame: returns the list of file paths that correspond to each group's fisher's exact test
    """
    feht_command="feht"
    try:
        if(shutil.which(snippy_command)) is not None:

            modified_data_frame.to_csv(output_directory+'/fishertest.txt', mode='w', index=False,sep="\t" )
            test_groups.to_csv(output_directory+'/testreference.txt', mode='w', index=False,sep="\t" )
            group_values=test_groups.group.unique()
            list_of_results=[]

            for i in range(len(group_values)):
            
                file_output=output_directory+"/feht_results"+group_values[i]+".txt"
                with open(file_output, 'w') as current_file:
                    argument_list=["feht", "-i", output_directory+"/testreference.txt", "-d", output_directory+"/fishertest.txt", "--one", "group "+ group_values[i], "-f", "1"]
                    call(argument_list, stdout=current_file)
                    list_of_results.append(file_output)
            
            return list_of_results 
        else:
            raise Exception("feht not installed on this machine, cannot proceed any further")
        
    except Exception as e:
        print(e.args)
        sys.exit(1)

    

   
    

    
