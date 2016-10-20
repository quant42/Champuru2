#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from simSeq import simulate as simSeqSimulate
from random import choice, randint, gauss

def choose(nuc, nucTerm):
    glob = nuc + nucTerm
    if glob == 0:
        return "E" # End
    r = randint(1, glob)
    if r > nuc:
        return "N" # Nuc
    return "T" # Terminal

def simPCRStrand(seq, prs):
    s = []
    for c in seq:
        nuc = choose(prs[c], prs[c + "*"])
        if nuc == "E":
            break
        elif nuc == "N":
            s.append(c)
            prs[c] = prs[c] - 1
        else:
            s.append(c)
            s.append("*")
            prs[c + "*"] = prs[c + "*"] - 1
            break
    return ("".join(s), prs)

def addSeqToResult(result, seq):
    result[seq] = result.get(seq, 0) + 1
    return result

def missPrs(prs):
    for key in prs:
        val = prs[key]
        if val <= 0:
            return True
    return False

def simulate():
    gcContent, len1, len2, startP1, startP2, snps, gap, testdata = simSeqSimulate()
    seqs = testdata[-1]
    assert len(seqs) == 4
    avogadro, dntp, ntp, acc = 6.022 * 10 ** 23, 0.5, 2, 0.01
    prs = {
        "A"  : int(round(gauss(ntp * avogadro, acc * avogadro))),
        "A*" : int(round(gauss(dntp * avogadro, acc * avogadro))),
        "C"  : int(round(gauss(ntp * avogadro, acc * avogadro))),
        "C*" : int(round(gauss(dntp * avogadro, acc * avogadro))),
        "T"  : int(round(gauss(ntp * avogadro, acc * avogadro))),
        "T*" : int(round(gauss(dntp * avogadro, acc * avogadro))),
        "G"  : int(round(gauss(ntp * avogadro, acc * avogadro))),
        "G*" : int(round(gauss(dntp * avogadro, acc * avogadro)))
    }
    resultFor = {}
    resultRev = {}
    for i in xrange(10000000):
        # not enough dntp's anymore? => quit
        if missPrs(prs):
            break
        # sim.
        i = randint(0, len(seqs) - 1) # the sequence index to take
        seq = seqs[i] # the sequence
        seq, prs = simPCRStrand(seq, prs)
        if i <= 1: # forward sequence
            resultFor = addSeqToResult(resultFor, seq)
        else: # reverse sequence
            resultRev = addSeqToResult(resultRev, seq)
    return (gcContent, len1, len2, startP1, startP2, snps, gap, testdata, resultFor, resultRev)

if __name__ == "__main__":
#    simulate()
    print(simulate())