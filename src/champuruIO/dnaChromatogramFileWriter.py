#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function

__doc__ = """
A dna chromatogram file writer.
"""

def writeToFile(filename, dic, comment):
    """
        Write chromatogram data to a chromatogram file.

        @param filename The name of the file to write the chromatogram to.
        @param dic      A dictionary giving the trace of each chromatogram.
        @param comment  A possible comment string to add to the chromatogram.
    """
    # check the file name. Until now only writing scf chromatogram files
    # are supported
    if not filename.lower().endswith(".scf"):
        filename = filename + ".scf"
    # ok, open that file
    with open(filename, "rb") as outF:
        pass
    pass # TODO
