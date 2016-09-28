#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import struct

__doc__ = """
A dna chromatogram file writer.
"""

def toBytes(l):
    result = b""
    for ele in l:
        if type(ele) == type(b""):
            result = result + ele
        elif type(ele) == type(0x00):
            result = result + chr(ele)
        else:
            assert False # should not go in here
    return result

def intToBytes(i, b): # little helper function
    if b == 4:
        # big-endian unsigned 4 encoding of int
        return struct.pack(">I", i)
    elif b == 2:
        # big-endian signed short encoding of int
        return struct.pack(">h", i)
    assert False # WTF? What am I expect to do?

def traceToByte(t):
    result, prev = [], 0
    for val in t:
        diff = val - prev
        result.append(intToBytes(diff, 2))
        prev = val
    return result

def writeToFile(filename, dic):
    """
        Write chromatogram data to a chromatogram file.

        @param filename The name of the file to write the chromatogram to.
        @param dic      A dictionary giving the trace of each chromatogram.
    """
    # check the file name. Until now only writing scf chromatogram files
    # are supported
    if not filename.lower().endswith(".scf"):
        filename = filename + ".scf"
    # check input length
    lengths = [len(dic[key]) for key in 'ACGT']
    for key in range(1, 4):
        assert lengths[0] == lengths[key], "Traces must be of same length"
    # ok, pre-calculate data
    ta = traceToByte(dic['A'])
    tc = traceToByte(dic['C'])
    tg = traceToByte(dic['G'])
    tt = traceToByte(dic['T'])
    assert len(ta) == len(tc) == len(tg) == len(tt), "Error while converting traces to bytes"
    # bases
    nrOfBases = 0
    # comment
    # means "PRG=Champuru"
    comment = [0x50, 0x52, 0x47, 0x3d, 0x43, 0x68, 0x61, 0x6d, 0x70, 0x75, 0x72, 0x75]
    # ok, open that file
    with open(filename, "wb") as outF:
        # write header
        # magic bytes (.scf)
        outF.write(toBytes([0x2e, 0x73, 0x63, 0x66]))
        # samples
        sampleLen = len(ta) // 2
        outF.write(intToBytes(sampleLen, 4))
        # sampleOffset
        outF.write(toBytes([0x00, 0x00, 0x00, 0x80]))
        # nrBases
        outF.write(intToBytes(len(dic['A']), 4))
        # nrBasesLeftClip
        outF.write(toBytes([0x00, 0x00, 0x00, 0x00]))
        # nrBasesRightClip - none
        outF.write(toBytes([0x00, 0x00, 0x00, 0x00]))
        # nrBasesOffset
        baseOffset = 128 + 4 * sampleLen
        outF.write(intToBytes(baseOffset, 4))
        # commentSize
        outF.write(intToBytes(len(comment), 4))
        # commentOffset
        commentOff = baseOffset + nrOfBases * 12
        outF.write(intToBytes(commentOff, 4))
        # version
        outF.write(toBytes([0x33, 0x2e, 0x30, 0x30]))
        # sampleSize
        outF.write(toBytes([0x00, 0x00, 0x00, 0x02]))
        # codeSet
        outF.write(toBytes([0x00, 0x00, 0x00, 0x00]))
        # privateSize - change if possible one day
        outF.write(toBytes([0x00, 0x00, 0x00, 0x00]))
        # privateOffset
        privOffset = commentOff + len(comment)
        outF.write(intToBytes(privOffset, 4))
        # spare
        outF.write(toBytes([0x00 for i in range(8 * 3)]))
        # write samples
        outF.write(toBytes(ta))
        outF.write(toBytes(tc))
        outF.write(toBytes(tg))
        outF.write(toBytes(tt))
        # write bases here - if possible one day
        # write comment
        outF.write(toBytes(comment))
        # write priv here if possible one day
