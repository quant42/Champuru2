#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import numpy as np
from scipy import signal

__doc__ = """
Library to convert the chromatogram into a probabilty matrix.
"""

def annotateMaxima(trace, traceCwt, pos):
    """
        Annotate a maxima.

        @param trace    The trace the maxima is in.
        @param traceCwt The continuous wavelet
        transformation of the trace.
        @param pos      The position of the maxima.
    """
    # get the value of the maximum
    maxVal = traceCwt[pos]
    # get the values of the nearest minimas
    lMinVal, lMinPos = maxVal, pos
    while lMinPos - 1 >= 0 and \
        lMinVal >= traceCwt[lMinPos - 1]:
        lMinPos -= 1
        lMinVal = traceCwt[lMinPos]
    rMinVal, rMinPos = maxVal, pos
    while rMinPos + 1 < len(traceCwt) and \
        rMinVal >= traceCwt[rMinPos + 1]:
        rMinPos += 1
        rMinVal = traceCwt[rMinPos]
    # calculate the quality of the maxima
    qual = 0 #TODO
    # return
    return (qual, pos)

def chromToMatrix(chrom):
    """
        Transfer the chromatogram into a probability
        matrix.

        @param chrom The chromatogram to process.
        @return A probability matrix that represents
        the amino acid probability for the different
        positions in the chromatogram.
    """
    # continuous wavelet transformation of chrom data
    cwt = {
        key : signal.cwt(chrom[key], signal.ricker, [4])[0]
        for key in "ACTG"
    }
    # Search for local maximums in the transformation
    maximas = {
        key : signal.argrelextrema(cwt[key], np.greater)[0]
        for key in "ACTG"
    }
    # annotate the maximas
    maximas = {
        key : [
            annotateMaxima(chrom[key], cwt[key], maxima)
            for maxima in maximas[key]
        ]
        for key in "ACTG"
    }
    # combine maximas
    
