import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def matread(filename):
      f = open(filename, 'r')
      l = [map(float,line.split()) for line in f]
      mat = np.array(l)
      return mat

c1 = matread("hydBbmat.csv")
c2 = matread("hydHbmat.csv")
c3 = matread("hydBHbmat.csv")

combined = np.concatenate((c1, c2, c3), axis=1)
#print combined.shape
#data = combined.transpose()
#print data.shape
name = ['L1','L2','L1-L2']


df = pd.DataFrame(data=combined, columns=name)
#print df

#df_norm = (df - df.mean())/(df.max() - df.min())

corr = df.corr(method='pearson')
print corr

#sns.pairplot(df)
#plt.show()

