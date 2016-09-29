#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from Bio._py3k import _bytes_to_string, _as_bytes
from dnaChromatogramFileReader import bytesToInt

__doc__ = """
Reader implementation for the scf file format specification.
"""

def isReadable(handle):
    """
        Fast check, if this file is readable by this reader. Check if the
        file magic bytes equals to '.scf' as specified in the file format
        specification:
        http://staden.sourceforge.net/manual/formats_unix_2.html

        @param handle The file handle.
        @return True if this (probably) is an abi file.
    """
    handle.seek(0)
    return handle.read(4) == _as_bytes('.scf')

def toSamples(data, sampleSize):
    # build the sample array
    samples = []
    while data != "":
        current, data = bytesToInt(data[0:sampleSize]), data[sampleSize:]
        if current & (1 << (sampleSize * 8 - 1)) != 0:   # check first bit (sign bit)
            current = -((1 << sampleSize * 8) - current) # convert to negative
        samples.append(current)
    # calculate delta
    for z in [1,2]:
        pSample = 0
        for i in range(len(samples)):
            samples[i] = samples[i] + pSample
            pSample = samples[i]
    return samples

def readFile(handle):
    """
        Read in the scf file.

        @param handle The file handle.
        @return The file data.
    """
    # go back to position 0
    handle.seek(0)
    # read in header
    header = {
        "magic" : handle.read(4),
        "samples" : bytesToInt(handle.read(4)),
        "sampleOffset" : bytesToInt(handle.read(4)),
        "nrBases" : bytesToInt(handle.read(4)),
        "nrBasesLeftClip" : bytesToInt(handle.read(4)),
        "nrBasesRightClip" : bytesToInt(handle.read(4)),
        "nrBasesOffset" : bytesToInt(handle.read(4)),
        "commentSize" : bytesToInt(handle.read(4)),
        "commentOffset" : bytesToInt(handle.read(4)),
        "version" : handle.read(4),
        "sampleSize" : bytesToInt(handle.read(4)),
        "codeSet" : handle.read(4),
        "privateSize" : bytesToInt(handle.read(4)),
        "privateOffset" : bytesToInt(handle.read(4)),
        "spare" : handle.read(4)
    }
    # check header
    if header["version"] != '3.00':
        raise Exception("Unsupported scf version! Until now only scf files with version 3.00 are supported. (Fileversion: `%s`)!" % header["version"])
    if header["sampleSize"] not in [1, 2]:
        raise Exception("Unknown scf sample size `%s`" % header["sampleSize"])
    # seek to sample offset
    handle.seek(header["sampleOffset"])
    # read in the traces => allways in order A, C, G, T for scf version 3.00
    aTrace, cTrace, gTrace, tTrace = \
        handle.read(header["samples"] * header["sampleSize"]), \
        handle.read(header["samples"] * header["sampleSize"]), \
        handle.read(header["samples"] * header["sampleSize"]), \
        handle.read(header["samples"] * header["sampleSize"])
    # check read in length
    if len(aTrace) != len(cTrace) != len(gTrace) != len(tTrace):
        raise Exception("Corrupted / Missing part in scf file `%s`!=`%s`!=`%s`!=`%s`" % (len(aTrace), len(cTrace), len(gTrace), len(tTrace)))
    # traces
    aTrace, cTrace, gTrace, tTrace = \
        toSamples(aTrace, header["sampleSize"]), \
        toSamples(cTrace, header["sampleSize"]), \
        toSamples(gTrace, header["sampleSize"]), \
        toSamples(tTrace, header["sampleSize"])
    # to interface
    data = {
        "FILE_EXT" : "SCF",
        "HEADER" : header,
        "DATA" : None, # TODO: not implemented
        "A" : aTrace,
        "C" : cTrace,
        "G" : gTrace,
        "T" : tTrace
    }
    # return the data
    return data
