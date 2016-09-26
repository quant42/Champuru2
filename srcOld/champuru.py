#! /usr/bin/env python
# -*- coding: UTF-8

from __future__ import division, print_function
import plotting, filter, cutter, processor, argparse
from champuruIO import dnaChromatogramFileReader as reader
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np

__doc__ = """
This file has no other purpose than to serve as main entry - in case you want to call
champuru from the commandline.
"""

if __name__ == "__main__":
    # argument parsing
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-i", "--forward", dest="in1", type=str, required=True, help="The chromatogram of the ''forward'' strand.")
#    parser.add_argument("-j", "--reversed", dest="in2", type=str, required=True, help="The chromatogram of the ''reversed'' strand.")
    args = parser.parse_args()
    # read in the chromatogram files
    print("Readin ...")
    i1 = reader.readFile(args.in1)
#    i2 = reader.readFile(args.in2)
#    plotting.plotChromatogram("i1.svg", i1)
#    plotting.plotChromatogram("i2.svg", i2)
    # filter
#    print("Filtering ...")
#    i1 = filter.applyGaussianFilter(i1)
#    i2 = filter.applyGaussianFilter(i2)
#    plotting.plotChromatogram("i1f.svg", i1)
#    plotting.plotChromatogram("i2f.svg", i2)
    print("Plotting ...")
    sign = signal.cwt(i1['A'], signal.ricker, [1, 50, 100, 500, 1000, 5000])
    for cwt in sign:
        plt.plot(cwt)
    plt.show()
    import sys; sys.exit(0)
    # cut out the reliable region out of the chromatogram
    i1, p11, p12 = cutter.cutoutReliableRegion(i1)
    i2, p21, p22 = cutter.cutoutReliableRegion(i2)
    plotting.plotChromatogram("i1c.svg", i1)
    plotting.plotChromatogram("i2c.svg", i2)
    # matrix representation
    m1 = processor.chromToMatrix(i1)
    m2 = processor.chromToMatrix(i2)
    print("Got sequences:")
    # TODO: represent matrices as seq.logo
    print(processor.matrixToSeq(m1))
    print(processor.matrixToSeq(m2))
    # reverse m2
    m2 = processor.reverse(m2)
    print("s2 reversed:")
    print(processor.matrixToSeq(m2))
    # combine matrix data
    offs, l1, l2 = [], len(m1), len(m2)
    minL = min(l1, l2)
    for offset in range(-l1+1, l2):
        minO, maxO = min(0, offset), max(0, offset)
        overlap = min(l1 + minO, l2 - maxO)
        mo1, mo2 = m1[-minO:-minO+overlap], m2[maxO:maxO+overlap]
        assert len(mo1) == len(mo2), "Internal error at (%s, %s)!" % (offset, overlap)
        score = processor.scoreOverlap(mo1, mo2)
        rScore = processor.getRandScore(mo1, mo2)
        rating = score - rScore
        offs.append((offset, score, rScore, overlap, rating))
    offs = sorted(offs, key=lambda x: x[4])
    print("Scores:")
    print(offs[:5])
    # TODO: plot offs stat
    # merge @ offset
    c1 = processor.combineMatrix(m1, m2, offs[0][0])
    c2 = processor.combineMatrix(m1, m2, offs[1][0])
    c3 = processor.combineMatrix(m2, m1, offs[0][0])
    c4 = processor.combineMatrix(m2, m1, offs[1][0])
    # reconstruct sequence
    print("Reconstructed sequences:")
    seq1 = processor.matrixToSeq(c1)
    seq2 = processor.matrixToSeq(c2)
    seq3 = processor.matrixToSeq(c1)
    seq4 = processor.matrixToSeq(c2)
    print(seq1)
    print(seq2)
    print(seq3)
    print(seq4)
