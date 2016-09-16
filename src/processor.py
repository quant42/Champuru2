#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import numpy as np
from scipy import signal, stats

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

def annotate(traceP, cwtP):
    """
        Annotate a part of a chromatogram part.

        @param traceP The part of the trace to
        annotate.
        @param cwtP   The part of the trace
        transformation to annotate.
        @return The annotation.
    """
    x = abs(max(0, max(cwtP)) - max(cwtP[0], cwtP[-1]))
    return 1 / (1 + np.e ** (-0.25*(x-20)))

def getPeakBetweenMinimas(chrom, start, stop):
    """
        Try to analyse/annotate a peak position.

        @param chrom The chromatogram the peak is in.
        @param start The start position in the
        chromatogram to search the peak in.
        @param stop  The stop position in the
        chromatogram to search the peak in.
        @return (start, stop, pA, pC, pT, pG) or
        None if between start and stop there's probably
        no peak!
    """
    # not enough data
    if stop - start <= 5:
        return None
    # ok, ok, analyse the data
    # TODO: check if tests like Kolmovorov-Smirnov-Test
    # Shapiro–Wilk work better
    data = []
    for key in "ACTG":
        data.append(annotate(chrom[key][start:stop+1], chrom['CWT_' + key][start:stop+1]))
    add = 1 - max(data)
    data = [data[i] + add for i in range(len(data))]
    return (start, stop, data[0], data[1], data[2], data[3])

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
    # calculate the continuous wavelet transformation
    # for each trace (for later)
    for key in "ACTG":
        chrom['CWT_' + key] = signal.cwt(chrom[key], signal.ricker, [4])[0]
    # Search for local minimas in the transformation
    minimas = signal.argrelextrema(mCwt, np.less)[0]
    # Now window between two consicutive minimas
    # this will explicitly not handle the case between
    # the last minima and the trace end
    startMin, lst = 0, []
    for minima in minimas:
        peak = getPeakBetweenMinimas(chrom, startMin, minima)
        if peak != None: # data got accepted
            lst.append(peak)
        startMin = minima
    # return the matrix
    return lst
