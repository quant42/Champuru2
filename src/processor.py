#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import numpy as np
from scipy import signal

__doc__ = """
Library to convert the chromatogram into a probabilty matrix.
"""

#def annotateMaxima(trace, traceCwt, pos):
#    """
#        Annotate a maxima.
#
#        @param trace    The trace the maxima is in.
#        @param traceCwt The continuous wavelet
#        transformation of the trace.
#        @param pos      The position of the maxima.
#        @return A triple (a, b, c), where a represents the
#        hight of the peak, b the quality of the peak and
#        c the peak position.
#    """
#    # get the value of the maximum
#    maxVal = traceCwt[pos]
#    # get the values of the nearest minimas
#    lMinVal, lMinPos = maxVal, pos
#    while lMinPos - 1 >= 0 and \
#        lMinVal >= traceCwt[lMinPos - 1]:
#        lMinPos -= 1
#        lMinVal = traceCwt[lMinPos]
#    rMinVal, rMinPos = maxVal, pos
#    while rMinPos + 1 < len(traceCwt) and \
#        rMinVal >= traceCwt[rMinPos + 1]:
#        rMinPos += 1
#        rMinVal = traceCwt[rMinPos]
#    # get the peak hight in the cwt transformation
#    cwtPeakHeight = maxVal - min(max(lMinVal, rMinVal), 0)
#    # calculate the quality of the maxima
#    # TODO better fit the parameter values of this model
#    qual = 1 / (1 + np.e ** (-0.25 * (cwtPeakHeight - 20)))
#    # return
#    return (cwtPeakHeight, qual, pos)

def chromToMatrix(chrom):
    """
        Transfer the chromatogram into a probability
        matrix.

        @param chrom The chromatogram to process.
        @return A probability matrix that represents
        the nucleotid probability for the different
        positions in the chromatogram.
    """
    # get the "maximal" trace
    maxTrace = []
    for i in range(len(chrom["A"])):
        maxVal = 0
        for key in "ACTG":
            maxVal = max(maxVal, chrom[key][i])
        maxTrace.append(maxVal)
    chrom['M'] = maxTrace
    # ok, in the maximal trace, search for local minimas
    # this allows us to reduce the problem, since between
    # each consicutive two local minima there's a local maxima
    # continuous wavelet transformation of this curve
    mCwt = signal.cwt(chrom['M'], signal.ricker, [4])[0]
    # Search for local minimas in the transformation
    minimas = signal.argrelextrema(mCwt, np.less)[0]
    # ok - now window between two consicutive minimas
    # this will explicitly not handle the case between
    # the last minima and the trace end
    startMin, lst = 0, []
    for minima in minimas:
        peak = getPeakBetweenMinimas(startMin, minima)
        lst.append(peak)
        startMin = minima
    # return the matrix
    return lst
