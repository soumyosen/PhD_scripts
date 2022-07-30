import numpy as np
import pandas as pd


number = []
with open("t.log") as f:
      for line in f:
            number.append(float(line))


#print(number)
#a=type(number)
#print(a)

data = np.reshape(number, (50,1))
#data1 = data.transpose()
#with open('hyd.BH.c1.mat.txt') as f:
#     for line in data1:
#         np.savetxt(f, line, fmt='%.2f')

dataframe = pd.DataFrame(data=data.astype(float))
dataframe.to_csv('hydHbmat.csv', sep=' ', header=False, float_format='%.2f', index=False)





#a = type(data1)
#print(data1)
#print(a)

#index = ['B4','B5','B6','B7','B8','B9']
#index = ['B10','B11','B12','B13','B14','B15','B16','B17']
#index = ['B18','B19','B20','B21','B22','B23','B24','B25']
#df=pd.DataFrame(data1, index=index)
#print(data1)
