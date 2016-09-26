#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from champuruIO import dnaChromatogramFileReader as reader
from champuruIO import dnaChromatogramFileWriter as writer
from scipy import signal
import numpy as np

__doc__ = """
This file basically consists out of the DNAChromatogram class
that is used in order to represent and handle DNA chromatograms.
"""

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
        # check if the lengths of each trace is equal
        # (or otherwise sayed, if there's any trace that has a
        # different length compared to the other traces).
        lengths = [len(aTrace), len(tTrace), len(cTrace), len(gTrace)]
        for i in range(1, len(lengths)):
            if lengths[0] != lengths[i]:
                raise Exception("Can't create new chromatogram object - traces differ in length")
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
    
    def save(self, filename, comment=None):
        """
            Save a chromatogram to a file.

            @param filename The filename to save this chromatogram to.
            @param comment  An additional string that should be saved
            as comment in the chromatogram file.
        """
        writer.writeToFile(
            filename,
            {
                'A' : self.aTrace,
                'C' : self.cTrace,
                'T' : self.tTrace,
                'G' : self.gTrace,
            },
            comment
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
    
    def plot(self, filename):
        """
            Get a graphical representation of this chromatogram.

            @param filename The filename to save the plot to.
        """
        pass # TODO
    
    def doBaseCalling(self):
        """
            Perform the base calling.

            @return 
        """
        # have we already done the base calling once?
        # if yes, just return the results
        try: return self.bases
        except: pass
        # nope, we haven't done the base calling till now
        # -> do it
        pass # TODO
    
    def getNucs(self):
        """
            Return the one letter code of the nucleotids that
            are treated with a DNA chromatogram. These can then
            be directly passed to chrom[key] in order to handle
            each trace independently.
        """
        return ['A', 'C', 'T' , 'G']
    
    def filter(self):
        """
            Filter the chromatogram.
        """
        pass # TODO
#TODO: annotate, __setitem__(self, key, val), iter(key, window=1)
