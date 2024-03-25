import pandas as pd
import numpy as np
folder = 'optimal_matrices/consumption/list_representation'
file1='TopologicalProperties_RandTrix_Nr25_Nc25_Nest0.4_Conn0.3312.csv'
file2='TopologicalProperties_RandTrix_Nr25_Nc25_Nest0.4_Conn0.3712.csv'

cols_to_use=[i for i in range(19) if i not in [7,8,14,16]]
data1 = pd.read_csv(folder+'/'+file1)
data2 = pd.read_csv(folder+'/'+file2)

string_keys = ['name', 'shared name']
keys = [a for a in pd.DataFrame.to_dict(data1).keys() if a not in string_keys]
strings = ['First column is for '+file1+'; second for '+file2]
for k in keys:
    string = k+' : '+str(data1[k].mean())+' +/- '+str(data1[k].std())+'; '+str(data2[k].mean())+' +/- '+str(data2[k].std())
    strings.append(string)
np.savetxt(folder+'/Analysis_TopologicalProperties.txt', strings,fmt="%s")
