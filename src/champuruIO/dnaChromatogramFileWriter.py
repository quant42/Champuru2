#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import struct

__doc__ = """
A dna chromatogram file writer.
"""

def probToByte(p):
    if p == None: return b"\x00"
    p = int(round(p * 255))
    assert 0 <= p < 256
    return intToBytes(p, 1)

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
    elif b == 1:
        return struct.pack("B", i)
    assert False # WTF? What am I expect to do?

def traceToByte(t):
    # convert to delta
    samples = t[:]
    for z in [1, 2]:
        pDelta = 0
        for i in range(len(samples)):
            pSample = samples[i]
            samples[i] = samples[i] - pDelta
            pDelta = pSample
    # convert deltas to bytes
    result = [intToBytes(val, 2) for val in samples]
    return result

def writeToFile(filename, dic, basecalling):
    """
        Write chromatogram data to a chromatogram file.

        @param filename    The name of the file to write the chromatogram to.
        @param dic         A dictionary giving the trace of each chromatogram.
        @param basecalling Basecalling data.
    """
    # check the file name. Until now only writing scf chromatogram files
    # are supported
    if not filename.lower().endswith(".scf"):
        filename = filename + ".scf"
    # check input length
    lengths = [len(dic[key]) for key in 'ACGT']
    for key in range(1, 4):
        assert lengths[0] == lengths[key], "Traces must be of same length"
    # convert basecalling input
    bases = []
    for b in basecalling:
        # ACTG
        nucMat = [ '-', 'A', 'C', 'M', 'T', 'W', 'Y', 'H', 'G', 'R', 'S', 'V', 'K', 'D', 'B', 'N' ]
        A, C, T, G = b[2:]
        index = (1 if A > 0.5 else 0) + \
            (2 if C > 0.5 else 0) + \
            (4 if T > 0.5 else 0) + \
            (8 if G > 0.5 else 0)
        bases.append(((b[0] + b[1]) // 2, b[2], b[3], b[5], b[4], nucMat[index]))
    assert len(bases) == len(set([x[0] for x in bases]))
    # ok, pre-calculate data
    ta = traceToByte(dic['A'])
    tc = traceToByte(dic['C'])
    tg = traceToByte(dic['G'])
    tt = traceToByte(dic['T'])
    assert len(ta) == len(tc) == len(tg) == len(tt), "Error while converting traces to bytes"
    # bases
    nrOfBases = len(bases)
    # comment
    # means "PRG=Champuru"
    comment = [0x50, 0x52, 0x47, 0x3d, 0x43, 0x68, 0x61, 0x6d, 0x70, 0x75, 0x72, 0x75]
    # ok, open that file
    with open(filename, "wb") as outF:
        # write header
        # magic bytes (.scf)
        outF.write(toBytes([0x2e, 0x73, 0x63, 0x66]))
        # samples
        sampleLen = len(ta) # ta = list of byte strings to write -> don't divide by 2
        outF.write(intToBytes(sampleLen, 4))
        # sampleOffset
        outF.write(toBytes([0x00, 0x00, 0x00, 0x80]))
        # nrBases
        outF.write(intToBytes(nrOfBases, 4))
        # nrBasesLeftClip
        outF.write(toBytes([0x00, 0x00, 0x00, 0x00]))
        # nrBasesRightClip - none
        outF.write(intToBytes(nrOfBases + 1, 4))
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
        outF.write(toBytes([0x00 for i in range(72)]))
        assert outF.tell() == 128
        # write samples
        outF.write(toBytes(ta))
        outF.write(toBytes(tc))
        outF.write(toBytes(tg))
        outF.write(toBytes(tt))
        assert outF.tell() == 128 + len(ta) * 2 * 4
        # write peak position
        for pos, a, c, t, g, s in bases:
            print("%s '%s'" % (pos, [str(hex(ord(b))) for b in intToBytes(pos, 4)]))
            outF.write(intToBytes(pos, 4))
        assert outF.tell() == 128 + len(ta) * 2 * 4 + nrOfBases * 4
        # acc.
        for i in xrange(1, 5):
            assert i in [1,2,3,4]
            for b in bases:
                outF.write(probToByte(b[i]))
        assert outF.tell() == 128 + len(ta) * 2 * 4 + nrOfBases * 8
        # sequence
        for base in bases:
            outF.write(base[5])
        assert outF.tell() == 128 + len(ta) * 2 * 4 + nrOfBases * 9
        # reserved
        outF.write(toBytes([0x00 for i in xrange(3 * nrOfBases)]))
        assert outF.tell() == 128 + len(ta) * 2 * 4 + nrOfBases * 12
        # write comment
        outF.write(toBytes(comment))
        # write priv here if possible one day
