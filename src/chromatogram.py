#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from champuruIO import dnaChromatogramFileReader as reader
from champuruIO import dnaChromatogramFileWriter as writer
from scipy.optimize import curve_fit
from scipy import signal
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
        if isinstance(aTrace, np.ndarray):
            aTrace = np.array(aTrace)
        if isinstance(cTrace, np.ndarray):
            cTrace = np.array(cTrace)
        if isinstance(tTrace, np.ndarray):
            tTrace = np.array(tTrace)
        if isinstance(gTrace, np.ndarray):
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
            return zip( # but programmed with numpy
                self.aTrace, self.cTrace, self.tTrace, self.gTrace
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
    
    def filter(self):
        """
            Filter the chromatogram.
        """
        def filterTrace(trace):
            # general
            fpoly = lambda x, a, b, c, d, e : a * x**4 + b * x**3 + c * x**2 + d * x + e
            # basic noise correction - low-pass filter
#            signal.butter(5, 100)
#            trace = signal.savgol_filter(trace, 5, 3)
            # get an estimation for the baseline
#            baseline, pcov = curve_fit(func, xdata, ydata)
            # get an estimation for the skyline
            cwtTrace = signal.cwt(trace, signal.ricker, [0.1])[0]
            maximas = signal.argrelextrema(cwtTrace, np.greater)[0]
            print(maximas)
            pSky, pCovSky = curve_fit(fpoly, maximas, np.array([trace[maxima] for maxima in maximas]))
            print(pSky)
            # give the filtered trace back
            return trace
        self.aTrace = filterTrace(self.aTrace)
        self.tTrace = filterTrace(self.tTrace)
        self.cTrace = filterTrace(self.cTrace)
        self.gTrace = filterTrace(self.gTrace)
    
    def getOverlapPosition(self, oChrom):
        """
            Get the overlap position, where the two chromatograms overlap.

            @param oChrom The chromatogram to overlap this chromatogram with. (This chromatogram
            should already be reversed, if it is representing a reversed chromatogram.)
            @return A list with indication positions.
        """
        def getTraceExtremas(chrom):
            result = {}
            for key in chrom.getNucs():
                # get the trace
                trace = chrom.__getitem__(key)
                # calculate the cwt transformation
                cwtTrace = signal.cwt(trace, signal.ricker, [0.1])[0]
                # get maximas in the cwt transformation
                maximas = signal.argrelextrema(cwtTrace, np.greater)[0]
                # save
                result[key] = maximas
            return result
        ext1 = getTraceExtremas(self)
        ext2 = getTraceExtremas(oTrace)
        # ok, try to best overlap the extremas
        
    
#TODO: annotate, __setitem__(self, key, val), iter(key, window=1)
# TODO: add base calling data to save in scf files
