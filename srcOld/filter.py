#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import numpy as np

__doc__ = """
Simply library for signal processing filters.
"""

def _getGaussianWeights(sigma):
    """
        Get the needed weights for the gaussian filter.

        @param sigma The sigma value for the gaussian
        filter to get the weights for.
        @return The needed weights for the gaussian
        filter.
    """
    weights = np.arange(0, 10*sigma+2, dtype=float) # after 10*sigma the weight is principally 0
    weights = (1 / (np.sqrt(2 * np.pi) * sigma)) * (np.e ** (-(weights ** 2 / (2 * sigma ** 2))))
    return weights

def gf(data, sigma):
    """
        Apply a gaussian filter to a signal.

        @param data  The signal to apply the filter on.
        @param sigma The sigma value of the gaussian filter.
        @return The filtered signal.
    """
    newData, weights, i = [], _getGaussianWeights(sigma), 0
    weightLen = min(len(weights), len(data))
    for pos, val in enumerate(data):
        summe = val * weights[0]
        for dist in range(1, weightLen):
            wd, pN, pP = weights[dist], pos - dist, pos + dist
            dpN = data[pN] if 0 <= pN < weightLen else 0
            dpP = data[pP] if 0 <= pP < weightLen else 0
            summe += dpN * wd
            summe += dpP * wd
        newData.append(int(round(summe)))
    return newData

def applyGaussianFilter(chrom):
    """
        Apply a gaussian filter on the chromatogram data.

        @param chrom The chromatogram to process.
        @return The processed chromatogram data.
    """
    # for all keys in the chromatogram
    for key in "ACTG":
        # get an estimation for the needed sigma value
        # Therefore use the end of the chromatogram where
        # the chromatogram normally isn't reliable anymore
        # in order to estimate the measurement fluctuation
        # error/noise of the sequencing machine ...
        # TODO: check if this can really be done everythime
        # e.g. expect that the user hasn't already cut of
        # the tail ...
        sigma = np.std(chrom[key][-100:], ddof=1)
        # apply filter for this trace
        chrom[key] = gf(chrom[key], sigma)
    # return the new chromatogram
    return chrom
