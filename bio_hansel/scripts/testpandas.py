import pandas as pd
data = [['Alex',10],['Bob',12],['Clarke',13]]
df = pd.DataFrame(data,columns=['Name','Age'])
new_data=df.drop('Name',1)
for index, row in df.iterrows():
        value_wanted=row['Name']
        print (df.iloc[1,1])

 print (df.iloc[1,1])