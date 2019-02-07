# -*- coding: utf-8 -*-

import math
import numpy as np 
import matplotlib.pyplot as plt 
import csv

# データ読み込み
p2 = np.genfromtxt('KarmanP.csv',delimiter=',')
p2 = p2.reshape([203,403])

# 描画
plt.contour(p2,50)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.show()
