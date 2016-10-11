#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from champuruIO import dnaChromatogramFileReader as reader
from champuruIO import dnaChromatogramFileWriter as writer
from scipy.optimize import curve_fit
from statistics import mean, stdev
from math import e, ceil, sqrt
from scipy import signal, fftpack
import numpy as np
import svgwrite
import random

__doc__ = """
This file basically consists out of the DNAChromatogram class
that is used in order to represent and handle DNA chromatograms.
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

def mms(seq):
    """
        Solver for the maximum sub sequence problem.

        @param seq The sequence to solve the maximum sub sequence for.
        @return Indices (start, stop) where the sequence is maximus.
    """
    # Algorithm - see: https://www.bio.ifi.lmu.de/mitarbeiter/volker-heun/notes/ab6.pdf
    maxscore, l, r = 0, 1, 0
    rmaxscore, rstart = 0, 1
    for i, val in enumerate(seq):
        if rmaxscore + val > val:
            rmaxscore += val
        else:
            rmaxscore, rstart = val, i
        if rmaxscore > maxscore:
            maxscore, l, r = rmaxscore, rstart, i
    return (l, r)

from itertools import islice

def window(seq, n=2): # see http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

class DNAChromatogram:
    """ Class representing a DNA chromatogram. """
    
    """ Measured values of the adenine trace """
    aTrace = np.array([])
    """ Measured values of the cytosine trace """
    cTrace = np.array([])
    """ Measured values of the thymine trace """
    tTrace = np.array([])
    """ Measured values of the guanine trace """
    gTrace = np.array([])
    """ The trace length """
    length = 0
    """ Whether the chromatogram got normalized """
    normalized = False
    
    def __init__(self, aTrace, cTrace, tTrace, gTrace):
        """
            Create a new chromatogram.

            @param aTrace The trace of adenine.
            @param cTrace The trace of cytosine.
            @param tTrace The trace of thymine.
            @param gTrace The trace of guanine.
        """
        # check the type of the traces; convert to np.array if needed
        # if this doesn't work -> error
        if not isinstance(aTrace, np.ndarray):
            aTrace = np.array(aTrace)
        if not isinstance(cTrace, np.ndarray):
            cTrace = np.array(cTrace)
        if not isinstance(tTrace, np.ndarray):
            tTrace = np.array(tTrace)
        if not isinstance(gTrace, np.ndarray):
            gTrace = np.array(gTrace)
        # check if there's no negative values in the traces
        assert min(aTrace) >= 0, "Negative measurement value in aTrace"
        assert min(cTrace) >= 0, "Negative measurement value in cTrace"
        assert min(tTrace) >= 0, "Negative measurement value in tTrace"
        assert min(gTrace) >= 0, "Negative measurement value in gTrace"
        # check if the lengths of each trace is equal
        # (or otherwise sayed, if there's any trace that has a
        # different length compared to the other traces).
        lengths = [len(aTrace), len(tTrace), len(cTrace), len(gTrace)]
        for i in range(1, len(lengths)):
            assert lengths[0] == lengths[i], "Can't create new chromatogram object - traces differ in length (%s)" % lengths
        self.length = lengths[0]
        # assign
        self.aTrace, self.cTrace, self.tTrace, self.gTrace = \
            aTrace, cTrace, tTrace, gTrace
    
    def reverse(self):
        """
            Reverse the chromatogram.
        """
        self.aTrace, self.cTrace, self.tTrace, self.gTrace = \
            self.tTrace[::-1], self.gTrace[::-1], self.aTrace[::-1], self.cTrace[::-1]
    
    @staticmethod
    def readFromFile(filename):
        """
            Create a DNA chromatogram object from a file.

            @param filename The file to read in.
            @return The chromatogram.
        """
        data = reader.readFile(filename)
        return DNAChromatogram(
            data['A'],
            data['C'],
            data['T'],
            data['G']
        )
    
    def save(self, filename):
        """
            Save a chromatogram to a file.

            @param filename The filename to save this chromatogram to.
        """
        writer.writeToFile(
            filename,
            {
                'A' : self.aTrace,
                'C' : self.cTrace,
                'T' : self.tTrace,
                'G' : self.gTrace,
            }
        )
    
    def __len__(self):
        """
            Return the length of the chromatogram traces.
        """
        return self.length
    
    def __getitem__(self, key):
        """
            Return the trace of a chromatogram by the trace key.

            @param key The trace to return the trace for.
        """
        if key == "A": return self.aTrace
        if key == "C": return self.cTrace
        if key == "T": return self.tTrace
        if key == "G": return self.gTrace
        if key == "Z":  # equivalent to the python zip transformation
            return np.array( # but programmed with numpy
                [self.aTrace, self.cTrace, self.tTrace, self.gTrace]
            ).T
        raise KeyError("Key: %s not defined for chromatogram object" % key)
    
    def __copy__(self):
        """
            Return a copy of the current object.
        """
        return DNAChromatogram(
            self.aTrace.__copy__(),
            self.cTrace.__copy__(),
            self.tTrace.__copy__(),
            self.gTrace.__copy__()
        )
    
    def getNucs(self):
        """
            Return the one letter code of the nucleotids that
            are treated with a DNA chromatogram. These can then
            be directly passed to chrom[key] in order to handle
            each trace independently.
        """
        return ['A', 'C', 'T' , 'G']
    
    def plot(self, filename, indents=0, colorCodes=[(255,0,0),(0,255,0),(0,0,255),(0,0,0)], keys='ACTG'):
        """
            Plot the chromatogram to a given file.

            @param filename The filename of the file to plot the
            chromatogram to.
            @param indents Number of indents values to plot at the beginning of the trace.
            @param colorCodes Color codes to use for plotting.
            @param keys The name/key of/for the traces to plot.
        """
        # check if filename ends with .svg
        if not filename.lower().endswith(".svg"):
            filename = filename + ".svg"
        # settings
        offsetX, offsetY = 50, 50
        # get the maximal value in the chromatogram
        maxChromVal = 0
        for key in keys:
            maxChromVal = max(maxChromVal, max(self.__getitem__(key)))
        # create an svg object
        svg = svgwrite.Drawing(filename=filename,
            size=(2 * offsetX + self.length + indents, 2 * offsetY + maxChromVal))
        # plot surrounding box
        drawbox = svg.rect(insert=(offsetX, offsetY), size=(self.length + indents, maxChromVal),
            fill="white", stroke="black")
        svg.add(drawbox)
        # plot legend
        for i, key in enumerate(keys):
            # add the box
            style = "fill:rgba(%s,%s,%s,0.3);" % colorCodes[i % len(colorCodes)]
            style += "stroke:rgb(%s,%s,%s)" % colorCodes[i % len(colorCodes)]
            box = svg.rect(insert=(offsetX + 15, offsetY + 15 + i * 24), size=(24, 12), style=style)
            svg.add(box)
            # add the text
            style = "font-size:20px"
            text = svg.text(key, insert=(offsetX + 52, offsetY + 27 + i * 24), style=style)
            svg.add(text)
        # plot bottom scale/ruler
        for i in range(0, self.length, 100):
            # ruler line
            line = svg.line(
                start=(offsetX + i + indents, offsetY + maxChromVal),
                end=(offsetX + i + indents, offsetY + maxChromVal + 10),
                style="stroke:black;stroke-width:1"
            )
            svg.add(line)
            # ruler text
            text = svg.text(str(i), insert=(offsetX + i + indents - 5, offsetY + maxChromVal + 28))
            svg.add(text)
        for i in range(0, int(ceil(maxChromVal)), 100):
            line = svg.line(
                start=(offsetX-10, offsetY + maxChromVal - i),
                end=(offsetX, offsetY + maxChromVal - i),
                style="stroke:black;stroke-width:1"
            )
            svg.add(line)
        # plot each trace
        for i, key in enumerate(keys):
            # plot the trace
            points = [(offsetX + indents, offsetY + maxChromVal)] # starting point is always "0", "0"
            lastPoint = offsetX + indents
            for j, value in enumerate(self.__getitem__(key)):
                lastPoint = offsetX + j + indents
                points.append((lastPoint, offsetY + maxChromVal - value))
            points.append((lastPoint, offsetY + maxChromVal)) # ending point is always "end", "0"
            # create the polygon
            style = "fill:rgba(%s,%s,%s,0.3);" % colorCodes[i % len(colorCodes)]
            style += "stroke:rgb(%s,%s,%s)" % colorCodes[i % len(colorCodes)]
            poly = svg.polygon(points, style=style)
            # add the polygon to the figure
            svg.add(poly)
        # save and return the svg object
        svg.save()
    
    def plotTraceAsColorMap(self, filename, traceKey, indents=0):
        """
            Plot the chromatogram to a given file. In contrast to the normal plot method,
            this plot only a trace, where the height data is represented by colors.
            This plot works best after baseline and skyline correction.

            @param filename The filename of the file to plot the
            chromatogram to.
            @param traceKey The key of the trace to plot.
            @param indents Number of indents values to plot at the beginning of the trace.
        """
        # check filename if it ends with .svg
        if not filename.lower().endswith(".svg"):
            filename = filename + ".svg"
        # settings
        offsetX, offsetY, height = 50, 50, 50
        # get the maximal value in the chromatogram
        maxChromVal = max(self.__getitem__(traceKey))
        # create an svg object
        svg = svgwrite.Drawing(filename=filename,
            size=(2 * offsetX + self.length + indents, 2 * offsetY + height))
        # plot surrounding box
        drawbox = svg.rect(insert=(offsetX, offsetY), size=(self.length + indents, height),
            fill="white", stroke="black")
        svg.add(drawbox)
        # plot bottom scale/ruler
        for i in range(0, self.length, 100):
            # ruler line
            line = svg.line(
                start=(offsetX + i + indents, offsetY + height),
                end=(offsetX + i + indents, offsetY + height + 10),
                style="stroke:black;stroke-width:1"
            )
            svg.add(line)
            # ruler text
            text = svg.text(str(i), insert=(offsetX + i + indents - 5, offsetY + height + 28))
            svg.add(text)
        # plot the trace
        trace = self.__getitem__(traceKey)
        for i, val in enumerate(trace):
            line = svg.line(
                start=(offsetX + i + indents, offsetY),
                end=(offsetX + i + indents, offsetY + height),
                style="stroke:hsl({},100%,50%);stroke-width:1".format(val / maxChromVal * 180)
            )
            svg.add(line)
        # save and return the svg object
        svg.save()
    
    def cutout(self, start, stop):
        """
            Get a chromatogram subpart.

            @param start The starting position.
            @param stop  The end position.
        """
        assert 0 <= start < stop <= self.length
        self.aTrace = self.aTrace[start:stop]
        self.cTrace = self.cTrace[start:stop]
        self.tTrace = self.tTrace[start:stop]
        self.gTrace = self.gTrace[start:stop]
        self.length = len(self.aTrace)
    
    def cutoutAuto(self, windowSize=100, treshold=20):
        """
            Cut out the uninterest start and end region from the chromatogram.

                - Combine the chromatogram traces into a single trace
                - Use a window of length [windowSize=100] to slide over the chromatogram.
                - Use fft to calculate the signal amplitude of each window position
                - Use mms to find where the signal is good.

            @param windowSize The size of the sliding over window.
            @param treshold A treshold (minimal amplitude) of the signal to analyse for filtering the rest of
            the trace.
        """
        # combine traces
        seq = [ sum(vals) for vals in self['Z'] ]
        # add a few values, so that seq has exactly the same length as the traces
        # this is needed in case, that the analysed sequence is longer than the chromatogram
        seq.extend([seq[len(seq) - 1 - i] for i in range(self.length + windowSize - 1 - len(seq))])
        assert len(seq) == self.length + windowSize - 1
        # use a window to iterate over the sequence (=combined trace)
        # and calculate the fft for each window
        amplis = []
        for pos, w in enumerate(window(seq, windowSize)):
            fft = fftpack.fft(w)
            amplitudes = 2.0 / 100 * np.abs(fft[:100//2])
            signalAmpli = max(amplitudes)
            amplis.append(signalAmpli)
        assert len(amplis) == self.length
#        import matplotlib.pyplot as plt
#        plt.plot(amplis)
#        plt.show()
        # convert to something useable by mms
        converged, start, stop = False, 0, len(amplis)
        while not converged:
            amplisF = filter(lambda x : x > treshold, amplis[start:stop])
            avg = mean(amplisF)
            std_ = 3 * stdev(amplisF)
            tSig = [-1 if abs(avg - val) > std_ or val < treshold else 1 for val in amplis]
            assert max(tSig) > 0
            # call mms problem solver
            startX, stopX = mms(tSig)
            startX = max(startX - windowSize // 2, 0)
            stopX = min(stopX + windowSize // 2, self.length)
            # check if converged
            converged = start == startX and stop == stopX
            start, stop = startX, stopX
        # cut
        self.cutout(start, stop)
        return (start, stop, amplis)
    
    def baseline(self):
        """
            Perform a baseline correction.
        """
        pass # basically this should already be done - TODO: check if this is really the case and implement an algorithm for threating this if this isn't the case
    
    def skyline(self, treshold=20):
        """
            Perform a skyline correction.

            @param treshold A treshold giving the minimal height of a peak to take into account for skyline correction.
        """
        def doSkylineCorrection(trace):
            # general
            tl = len(trace)
            # skyline estimation
            cwtTrace = signal.cwt(trace, signal.ricker, [0.1])[0]
            maximas = signal.argrelextrema(cwtTrace, np.greater)[0]
            maxVals = [(maxima, cwtTrace[maxima]) for maxima in maximas] # this is for maxima filtering
            maxVals = filter(lambda x:x[1]>treshold, maxVals) # TODO: check if there's a better treshold
            # create the maxima arrays
            maximas = np.array([x[0] for x in maxVals])
            maxVals = np.array([x[1] for x in maxVals])
            # expected peak height ~ exp. decay ratio
            pSkyF = lambda x, a, b : a * e ** (-b * x)
            (a, b), pConv = curve_fit(pSkyF, maximas, maxVals, p0=(max(maxVals), 10E-7))
            skyline = np.array([pSkyF(x, a, b) for x in range(tl)])
            skyline = np.array([max(1, skyline[x]) for x in range(tl)])
            # normalize
            trace = np.array([ trace[x] / skyline[x] * 100 for x in range(tl) ])
            # return
            return trace
        self.aTrace = doSkylineCorrection(self.aTrace)
        self.cTrace = doSkylineCorrection(self.cTrace)
        self.tTrace = doSkylineCorrection(self.tTrace)
        self.gTrace = doSkylineCorrection(self.gTrace)
    
    def noiseCorrection(self):
        """
            Perform a noise correction (filtering).
        """
        pass # Same as with baseline - This should already be done. TODO: check if this is really the case and do it if not.
    
    def doBaseCalling(self, params=(1.61, 0.1, 6, 1.38, 12)):
        """
            Perform a base calling. (Only do that after cutoutAuto, baseline, skyline and noiseCorrection have
            been applied.)

            @param params Basecalling parameters.
        """
        # get the "maximal" trace
        maxTrace = []
        for s in self['Z']:
            maxVal = max(s)
            maxTrace.append(maxVal)
        # ok, in the maximal trace, search for local minimas
        # this allows us to reduce the problem, since between
        # each consicutive two local minima there's a local maxima
        # continuous wavelet transformation of this curve
        mCwt = signal.cwt(maxTrace, signal.ricker, [params[0]])[0]
        # calculate the continuous wavelet transformation
        # for each trace (for later)
        for key in self.getNucs():
            cwt = signal.cwt(self[key], signal.ricker, [params[1]])[0]
            setattr(self, "CWT_%s" % key, cwt)
        # Search for local minimas in the transformation
        minimas = signal.argrelextrema(mCwt, np.less)[0]
        # assure that the end point of trace is in
        if minimas[-1] < self.length - 2:
            minimas = np.append(minimas, self.length - 1)
        # some helper functions
        def annotate(traceP, cwtP, params):
            x = abs(max(0, max(cwtP)) - max(cwtP[0], cwtP[-1]))
            # prevent overflow
            if (-params[3]*(x-params[4])) < -100: return 1.0
            if (-params[3]*(x-params[4])) > 100: return 0
            # calculate precisely
            return 1 / (1 + np.e ** (-params[3]*(x-params[4])))
        # Now window between two consicutive minimas
        startMin, lst = 0, []
        for minima in minimas:
            # convinience renaming
            start, stop = startMin, minima
            # check if enough data
            if stop - start <= params[2]:
                lst.append((start, stop, None, None, None, None))
            peakData = []
            for key in self.getNucs():
                peakData.append(
                    annotate(self[key][start:stop+1], getattr(self, 'CWT_%s' % key)[start:stop+1], params)
                )
            lst.append((start, stop, peakData[0], peakData[1], peakData[2], peakData[3]))
            startMin = minima
        # return the matrix
        return lst
    
    def peakIdentifing(self):
        """
            Identify peaks in the chromatogram for the different traces.

            @return An dictionary ACTG with the corresponding peak position.
        """
        result = {}
        for key in self.getNucs():
            mCwt = signal.cwt(self[key], signal.ricker, [0.1])[0]
            maximas = signal.argrelextrema(mCwt, np.greater)[0]
            result[key] = maximas
        return result
    
    def getFuzzyRep(self):
        """
            Get a fuzzy representation of this chromatogram.
        """
        if not self.normalized: self.normalize()
        def valToFuzzy(val): # return value between inclusive 0 and 5
            assert val >= 0
            return min(val // 20, 5)
        seq = []
        for va, vc, vt, vg in self['Z']:
            seq.append((valToFuzzy(va), valToFuzzy(vc), valToFuzzy(vt), valToFuzzy(vg)))
        return np.array(seq)
    
    def align(self, oChrom):
        """
            Align this chromarogram with another chromatogram.

            @param oChrom the other chromatogram to align this chromatogram with.
        """
        # get a string representation of the chromatogram
        s1 = self.getFuzzyRep()
        s2 = oChrom.getFuzzyRep()
        # create the matrix
        matrix = np.zeros((len(s1) + 1, len(s2) + 1))
        # fill out the matrix
        def d(a, b): # val between inclusive -2 and 3
            return sqrt(a * b) - 2
        def diff(v1, v2):
            return max([d(v1[i], v2[i]) for i in range(4)])
        for p1, v1 in enumerate(s1):
            for p2, v2 in enumerate(s2):
                matrix[p1+1][p2+1] = max(
                    0, # local alignment
                    matrix[p1][p2+1] - 1,
                    matrix[p1+1][p2] - 1,
                    matrix[p1][p2] + diff(v1, v2)
                )
        # TODO: heatmap of matrix
        import matplotlib.pyplot as plt
        plt.pcolor(matrix)
        plt.savefig("matrix.png")
        # find maximas in matrix
        
#TODO: annotate, __setitem__(self, key, val), iter(key, window=1)
# TODO: add base calling data to save in scf files
