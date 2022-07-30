import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def matread(filename):
      f = open(filename, 'r')
      l = [map(float,line.split()) for line in f]
      mat = np.array(l)
      return mat

c1 = matread("hydBHc1mat.csv")
c = matread("hydBHcmat.csv")
c2 = matread("hydBHc2mat.csv")

combined = np.concatenate((c1, c, c2), axis=0)
#print combined.shape
data = combined.transpose()
#print data.shape
name = ['B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25']


df = pd.DataFrame(data=data, columns=name)
#print df

corr = df.corr(method='pearson')
#print corr

sns.pairplot(df)
plt.show()

