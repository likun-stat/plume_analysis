#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 22:10:45 2024

@author: LikunZhang
"""

import os

import csv
with open('./xz_all.csv', newline='') as csvfile:
    X_csv = list(csv.reader(csvfile, delimiter=' ', quotechar='|'))


import numpy as np
X = np.empty((len(X_csv)-1, 900))
for iter in np.arange(1,len(X_csv)):
    tmp = X_csv[iter][0].split(",")
    X[iter-1,:] = [float(i) for i in tmp[1:901]]


from sklearn.decomposition import NMF
# 7 13 19 25 31 37 43 48 54 60 66 72 78 84 90
n_comp = 50
model = NMF(n_components=n_comp, init='random', random_state=0, max_iter=10000)
W = model.fit_transform(X)
H = model.components_


# import seaborn as sns
# import matplotlib.pylab as plt

# w_data = W[:,1].reshape(198,500)
# ax = sns.heatmap(w_data, linewidth=0.5)
# plt.show()

np.savetxt("features"+str(n_comp)+".csv", W, delimiter = ",")
np.savetxt("components"+str(n_comp)+".csv", H, delimiter = ",")
