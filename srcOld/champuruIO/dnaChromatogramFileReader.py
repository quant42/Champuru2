#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import sys, os, importlib

__doc__ = """
A reader, that is allowing to read in the different DNA chromatogram file
formats.
"""

# The files containing dna chromatogram reader implementation for different
# file types each reader needs to implement the two methods isReadable and
# readFile.
global implementedReaders
implementedReaders = [
    "abiReader.py", # reader for abi format
    "scfReader.py", # reader for scf format
]

def readFile(filename):
    """
        Read the data from an DNA chromatogram file. This method is checking
        if there's a reader implementation that allows to read in this DNA
        chromatogram data. If such a reader implementation get's found, this
        reader implementation will get used in order to read in the specific
        file.

        @param filename The filename of the file to read in.
        @return A dictionary representing the data of the file.
        @throws Exception, if there's no reader that can read this file, or
        if something goes wrong while reading in the file data.
    """
    # open the file - only open the file once. This is faster.
    with open(filename, "rb") as inF:
        # check if there's a reader implemation, that is able to read the file
        for reader in implementedReaders:
            # seek to 0. This is operation may be needed, if the the previous
            # readers isReadable method messed up the file position.
            # Next reader may expect that reading starts at file position 0.
            inF.seek(0)
            # import the reader module - out of the current directory of course
            sys.path.insert(0, os.path.dirname(__file__))
            try:
                readerImpl = importlib.import_module(reader[:-3])
            finally:
                sys.path.pop(0)
            # check if the reader chan handle this file
            if readerImpl.isReadable(inF):
                # cool, probably this reader can handle this file
                # so let this reader handle this file
                inF.seek(0)
                return readerImpl.readFile(inF)
    # no reader supported this file -> raise exception
    raise Exception("Unsupported file format - Can't read data of file '{}'."
        .format(filename))

# functions that can be used by all readers
bytesToInt = lambda x : int(str(x).encode('hex'), 16)
