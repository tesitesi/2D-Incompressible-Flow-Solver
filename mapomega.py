# -*- coding: utf-8 -*-

import math
import numpy as np 
import matplotlib.pyplot as plt 
import csv

# データ読み込み
p2 = np.loadtxt('KarmanOmega.csv',delimiter=',')

# データ２次元配列作成

# 描画
plt.contour(p2,50)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.show()
