#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from chromatogram import DNAChromatogram as chromo

d = chromo.readFromFile("/home/quant/data/binf/chromatograms/abi/ExampleForward.ab1")
#e = chromo.readFromFile("/home/quant/data/binf/chromatograms/abi/ExampleReverse.ab1")

d.cutoutAuto()

#e.cutoutAuto()
#d.noiseCorrection()
#e.noiseCorrection()
#d.baseline()
#e.baseline()
#d.skyline()
#e.skyline()

#e.reverse()

d.plot("frw.svg")
#e.plot("rev.svg")
import sys
sys.exit(0)
dPeaks = d.peakIdentifing()
ePeaks = e.peakIdentifing()

def getDists(l):
    for i, ele in enumerate(l):
        if i == 0: continue
        yield ele - l[i-1]

dPeakDist = {k:list(getDists(l)) for k,l in dPeaks.items()}
ePeakDist = {k:list(getDists(l)) for k,l in ePeaks.items()}
print(dPeaks)
print(dPeakDist)
