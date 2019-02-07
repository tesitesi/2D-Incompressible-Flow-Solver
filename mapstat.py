# -*- coding: utf-8 -*-

import math
import numpy as np 
import matplotlib.pyplot as plt 
import csv

# データ読み込み
data = np.loadtxt('KarmanStat.csv',delimiter=',')

n = data[:,0]
res = data[:,1]
cd = data[:,2]
cl = data[:,3]
cp1 = data[:,4]
cp2 = data[:,5]

plt.plot(n,cd,label="CD")
plt.plot(n,cl,label="CL")
plt.plot(n,cp1,label="Cp1")
plt.plot(n,cp2,label="Cp2")
plt.xlabel('No. of Steps')
plt.legend()
plt.show()
