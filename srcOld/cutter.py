#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function

__doc__ = """
A simple library for cutting of the non reliable beginning
and end of a filtered chromatogram.
"""

def getRCutPos(chrom):
    """
        Get the right position where the chromatogram
        should get cut of.

        @param chrom The chromatogram to process.
        @return The positon.
    """
    # the position to cut of. (-1 indicates that no position
    # got found)
    rCutPos = -1
    # see which trace get's up last
    for key in "ACTG":
        pos = len(chrom[key]) - 1
        # after filter, the noise level should be reduced
        # to a minimum ...
        while pos > 0 and chrom[key][pos] <= 4:
            pos -= 1
        rCutPos = max(rCutPos, pos)
    # return the position or the chromatogram length if no
    # such positon has been found.
    return rCutPos if rCutPos != -1 else len(chrom["A"])

def getLCutPos(chrom):
    """
        Get the left position where the chromatogram
        should get cut of.

        @param chrom The chromatogram to process.
        @return The positon.
    """
    lCutPos, flag = 0, False
    # curves should be reliable at N//2 ->
    # so anyway only check till that position
    for pos in range(2, len(chrom["A"]) // 2):
        count = 0
        for key in "ATCG":
            summe, falling = 0, True
            for okey in "ATCG":
                if key != okey:
                    summe += chrom[okey][pos]
                    arr = chrom[okey][pos:pos+15]
                    if max(arr) != arr[0]:
                        falling = False
                        break
            if falling and \
                1.5 * summe < chrom[key][pos] and \
                6 < chrom[key][pos]:
                lCutPos = pos
        if lCutPos != 0: break
    return lCutPos

def cutoutReliableRegion(chrom):
    """
        Cut off the unreliable chromatogram region at
        the beginning and the end of a *filtered*
        chromatogram.

        @param chrom The chromatogram to process.
        @return A triple (c, p1, p2) where c is the
        trimmed chromatogram and p1, p2 are representing
        the cutting positions.
    """
    # get the cutting positions
    lPos, rPos = getLCutPos(chrom), getRCutPos(chrom)
    # chop
    for key in "ATCG":
        chrom[key] = chrom[key][lPos:rPos]
    # return
    return (chrom, lPos, rPos)
