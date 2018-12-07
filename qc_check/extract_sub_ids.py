'''
This script reads the csv file of paths of bugs and returns the array of bugs.
'''
import pandas as pd

in_file = 'ABIDE_2_Bug_paths.csv'

df = pd.read_csv(in_file).values

count = 0
bugs = []
for idx in range(df.shape[0]):
    sub_id = df[idx,0].split('/')[-1].split('_')[0].split('-')[1]
    # print(sub_id)
    count = count + 1
    bugs.append(sub_id)

print(bugs)
print(count)
