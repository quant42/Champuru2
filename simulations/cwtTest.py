#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from scipy.optimize import curve_fit
import numpy as np
from scipy import signal, fftpack
from math import sqrt, e, pi
import matplotlib.pyplot as plt
from itertools import islice
from statistics import mean as avg
from statistics import stdev as std

def window(seq, n=2): # see http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-py$
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def findMaximas(trace):
    maximas = []
    for pos, vals in enumerate(window(trace, 5)):
        maxEle = max(vals)
        if vals[2] == maxEle and vals[1] < vals[2] >= vals[3]:
            maximas.append(pos + 2)
    return maximas

gauss = lambda x, mu, sigma : 1 / (sqrt(2 * pi * sigma ** 2)) * e ** (- (x - mu) ** 2 / (2 * sigma ** 2))
ff = lambda x, a1, mu1, sigma1, a2, mu2, sigma2 : a1 * gauss(x, mu1, sigma1) + a2 * gauss(x, mu2, sigma2)

n, sigma = 100, 0.2
coord = np.linspace(4, 7, n)
#print(list(coord))
f = [ff(x, 20, 5, sigma, 6, 5.5, sigma) for x in coord]
f_ = [20 * gauss(x, 5, sigma)  for x in coord]
f__ = [6 * gauss(x, 5.5, sigma)  for x in coord]
plt.plot(coord, f, label="orig")
plt.plot(coord, f_, label="orig")
plt.plot(coord, f__, label="orig")

#print(curve_fit(ff, coord, f)[0])
print([coord[x] for x in findMaximas(f)])
#print(curve_fit(gauss, coord, f)[0])

ws = [2.5*2] #[10, 5, 0.1]
for k, mwlet in enumerate([signal.ricker]): #, signal.gaussian]):
    cwts = signal.cwt(f, mwlet, ws)
    for i, cwt in enumerate(cwts):
        plt.plot(coord, cwt, label="cwt_%s_%s" % (k, ws[i]))
        print([coord[x] for x in findMaximas(cwt)])
plt.legend(loc="best")

plt.show()
