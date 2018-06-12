import pandas as pd
from subprocess import call


def conductFisherTest(data_frame):
    with open('print.sh', 'rb') as file:
      script = file.read()
    modified_df=data_frame
    modified_df.to_csv('fishertest.txt', mode='w', index=False,sep="\t" )

    # rc = call(script, shell=True)
    

    
