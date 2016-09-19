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

def matrixToSeq(m):
    """
        Transfer a matrix into a sequence.
    """
    result = []
    for data in m:
        add = 1 - max(data[:4])
        A = data[0] + add > 0.5
        C = data[1] + add > 0.5
        T = data[2] + add > 0.5
        G = data[3] + add > 0.5
        val = (1 if A else 0) + \
            (2 if C else 0) + \
            (4 if T else 0) + \
            (8 if G else 0)
        result.append([
            'Z', 'A', 'C', 'M', 'T', 'W', 'Y', 'H', 'G',
            'R', 'S', 'V', 'K', 'D', 'B', 'N'
        ][val])
    return "".join(result)

def combine(c1, c2, p1, p2):
    """
        Combine two positions in the matrix.

        @param c1 The first matrix.
        @param c2 The second matrix.
        @param p1 The position in the first matrix.
        @param p2 The position in the second matrix.
        @result A combined probability vector.
    """
    try:
        return ( # ACTGXY
            (c1[p1][0] * c2[p2][2]) ** 0.5,
            (c1[p1][1] * c2[p2][3]) ** 0.5,
            (c1[p1][2] * c2[p2][0]) ** 0.5,
            (c1[p1][3] * c2[p2][1]) ** 0.5
        )
    except:
        if 0 <= p1 < len(c1):
            return c1[p1]
        if 0 <= p2 < len(c2):
            return c2[p2]
        return (0,0,0,0)

def combineMatrix(c1, c2, offset):
    """
        Combine two matrices.

        @param c1 The first matrix.
        @param c2 The second matrix.
        @param offset The offset for combination.
        @return The combined matrix.
    """
    return [
        combine(c1, c2, i - offset, i)
        for i in range(min(offset, 0), max(len(c1), len(c2)) + max(offset, 0))
    ]

def evaluatePairing(c1, c2):
    """
        Evaluate a pairing.

        @param c1 The first of the pairings.
        @param c2 The second of the pairings.
        @result A pairing score between 0 (poor pairing)
        and 1 (very good pairing).
    """
    val = max(
        (c1[0] * c2[2]) ** 0.5,
        (c1[1] * c2[3]) ** 0.5,
        (c1[2] * c2[0]) ** 0.5,
        (c1[3] * c2[1]) ** 0.5
    )
    assert 0 <= val <= 1, "Internal Error; 0 <= %s <= 1?" % val
    return val

def scoreOverlap(c1, c2):
    """
        Score an overlap.

        @param c1 The overlap of the first chromatogram.
        @param c2 The overlap of the second chromatogram.
        @return The score of the overlaps.
    """
    score = 0
    for i, cS1 in enumerate(c1):
        score += evaluatePairing(cS1, c2[i])
    return score

def getRandScore(c1, c2):
    """
        Get the score, that a random overlap, with the base
        distribution of c1 and c2 would probably archive.

        @param c1 The overlap of the first chromatogram.
        @param c2 The overlap of the second chromatogram.
        @return The score that a random overlap would archive.
    """
    summe = 0
    for ch1 in c1:
        for ch2 in c2:
            summe += evaluatePairing(ch1, ch2)
    return summe / len(c1)

def reverse(matrix):
    """
        Reverse a matrix.

        @param matrix The matrix to return.
        @return The reversed matrix.
    """
    mRet = [ # ACTG
        (val[2], val[3], val[0], val[1])
        for val in reversed(matrix)
    ]
    return mRet

def annotate(traceP, cwtP, params):
    """
        Annotate a part of a chromatogram part.

        @param traceP The part of the trace to
        annotate.
        @param cwtP   The part of the trace
        transformation to annotate.
        @param params Annotation parameters.
        @return The annotation.
    """
    x = abs(max(0, max(cwtP)) - max(cwtP[0], cwtP[-1]))
    # prevent overflow
    if (-params[3]*(x-params[4])) < -100: return 1.0
    if (-params[3]*(x-params[4])) > 100: return 0
    # calculate precisely
    return 1 / (1 + np.e ** (-params[3]*(x-params[4])))

def getPeakBetweenMinimas(chrom, start, stop, params):
    """
        Try to analyse/annotate a peak position.

        @param chrom The chromatogram the peak is in.
        @param start The start position in the
        chromatogram to search the peak in.
        @param stop  The stop position in the
        chromatogram to search the peak in.
        @param params Prediciton parameters.
        @return (pA, pC, pT, pG, start, stop) or
        None if between start and stop there's probably
        no peak!
    """
    # not enough data
    if stop - start <= params[2]:
        return None
    # ok, ok, analyse the data
    # TODO: check if tests like Kolmovorov-Smirnov-Test
    # Shapiro–Wilk work better
    data = []
    for key in "ACTG":
        data.append(annotate(chrom[key][start:stop+1], chrom['CWT_' + key][start:stop+1], params))
#    add = 1 - max(data)
#    data = [data[i] + add for i in range(len(data))]
    return (data[0], data[1], data[2], data[3], start, stop)

def chromToMatrix(chrom, params=(1.61, 0.1, 6, 1.38, 12)):
    """
        Transfer the chromatogram into a probability
        matrix.

        @param chrom The chromatogram to process.
        @param params Possible prediction parameters.
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
    mCwt = signal.cwt(chrom['M'], signal.ricker, [params[0]])[0]
    # calculate the continuous wavelet transformation
    # for each trace (for later)
    for key in "ACTG":
        chrom['CWT_' + key] = signal.cwt(chrom[key], signal.ricker, [params[1]])[0]
    # Search for local minimas in the transformation
    minimas = signal.argrelextrema(mCwt, np.less)[0]
    # Now window between two consicutive minimas
    # this will explicitly not handle the case between
    # the last minima and the trace end
    startMin, lst = 0, []
    for minima in minimas:
        peak = getPeakBetweenMinimas(chrom, startMin, minima, params)
        if peak != None: # data got accepted
            lst.append(peak)
        startMin = minima
    # return the matrix
    return lst
