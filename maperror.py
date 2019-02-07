# -*- coding: utf-8 -*-

import math
import numpy as np 
import matplotlib.pyplot as plt 
import csv

# データ読み込み
data = np.loadtxt('KarmanStat.csv',delimiter=',')

n = data[:,0]
res = data[:,1]


plt.plot(n,res)
plt.xlabel('No. of Steps')
plt.ylabel('RMS Error for p')
plt.yscale("log")
plt.show()
