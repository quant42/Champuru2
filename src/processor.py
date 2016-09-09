#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import statistics
import numpy as np
from scipy import signal

__doc__ = """
Library to convert the chromatogram into a probabilty matrix.
"""

def getHighlyReliablePeaks(chrom):
    """
        Get high reliable peaks in a chromatogram.

        @chrom The chromatogram to analyse.
        @return The positions (sorted) of the highly reliable
        peaks in the chromatogram.
    """
    positions = set()
    for key in "ACTG":
        # TODO: find_peaks_cwt is not optimal for this case ...
        # find or implement a better method (or add a transformation
        # that is making peaks_cwt better suited)
        posis = signal.find_peaks_cwt(chrom[key], np.arange(1,10))
        positions.update(posis)
    return sorted(list(positions))

def getPeakDistance(chrom):
    """
        Get the average distance between two peaks.

        @param chrom The chromatogram to analyse.
        @return The average distance between two
        peaks.
    """
    # get the positions of the hightly reliable peaks
    # in the chromatogram
    peaks = getHighlyReliablePeaks(chrom)
    # calculate the distance between the peaks
    peakDists = []
    for i, peak in enumerate(peaks):
        if i != 0:
            peakDists.append(peak - peaks[i-1])
    # It is expectable that the meadian distance between
    # the reliable peaks is the average peak distance
    return statistics.median(peakDists)

def chromToMatrix(chrom):
    """
        Transfer the chromatogram into a probability
        matrix.

        @param chrom The chromatogram to process.
        @return A probability matrix that represents
        the amino acid probability for the different
        positions in the chromatogram.
    """
    # first step: get the distance between the peaks
    dist = getPeakDistance(chrom)
