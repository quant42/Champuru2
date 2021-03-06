#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from simSeq import simulate as simSeqSimulate
from random import choice, randint, gauss, getrandbits, random, randrange
try: xrange
except NameError: xrange = range

def doPCR(seq, prs):
    flag, s = False, []
    for c in seq:
        cS = c + "*"
        a, b = prs[c], prs[cS]
        glob = a + b
        if glob == 0:
            flag = True
            break
        r = randrange(glob)
        if r < a:
            s.append(cS)
            prs[cS] = prs[cS] - 1
            break
        s.append(c)
        prs[c] = prs[c] - 1
    return (flag, "".join(s))

def simulate():
    gcContent, len1, len2, startP1, startP2, snps, gap, testdata = simSeqSimulate()
    seqs = testdata[-1]
    assert len(seqs) == 4
    avogadro, dntp, ntp, acc = 6.022 * 10 ** 23, 1/(2 ** len(seqs[0])), 1, 0.01
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
    for i in xrange(1000000):
        # sim.
        i = getrandbits(2) # fast way for writing randint(0, 3) - this function is not cryptographically secure, but that doesn't mind for PCR simulations. # the sequence index to take; 2 * 2 (diploid * (forw+rev))
        seq = seqs[i] # the sequence
        flag, seq = doPCR(seq, prs)
        # not enough dntp's anymore? => quit
        if flag:
            break
        if i <= 1: # forward sequence
            resultFor[seq] = resultFor.get(seq, 0) + 1
        else: # reverse sequence
            resultRev[seq] = resultRev.get(seq, 0) + 1
    return (gcContent, len1, len2, startP1, startP2, snps, gap, testdata, resultFor, resultRev)

if __name__ == "__main__":
#    simulate()
    print(simulate())
