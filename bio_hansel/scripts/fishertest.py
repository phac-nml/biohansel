import pandas as pd
from subprocess import call


def conductFisherTest(modified_data_frame, output_directory, test_groups):

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
    

    # rc = call(script, shell=True)
    

    
