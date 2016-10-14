#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function

"""
Some general helper stuff about nucleotids.
"""

# DATA
_nc2code = [ # map a code to each character
#  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15
 '?', 'A', 'T', 'W', 'G', 'R', 'K', 'D', 'C', 'M', 'Y', 'H', 'S', 'V', 'B', 'N'
]
_pairs = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C' } # the nucleotids that "pairs"

# haskel like helper functions
def splitBin(s):
    result, i = [], 0
    while s > 0:
        x = 1 << i
        if s & x != 0:
            result.append(x)
            s = s ^ x
        i += 1
    return result

charToCode = lambda c : _nc2code.index(c)
strToCodes = lambda s : [charToCode(c) for c in s]
charsToCode = lambda l : sum([charToCode(c) for c in l])
codeToChar = lambda c : _nc2code[c]
codeToChars = lambda c : [codeToChar(s) for s in splitBin(c)]
codeToECode = lambda c : charsToCode([_pairs[codeToChar(x)] for x in splitBin(c)])
charsToChar = lambda c : codeToChar(charsToCode(c))
codesToStr = lambda l : "".join([codeToChar(c) for c in l])
reverseStr = lambda s : "".join([codeToChar(codeToECode(charToCode(c))) for c in s[::-1]])

# filled out fast lookups
nc2code = {j : i for i, j in enumerate(_nc2code)}
code2nc = {i : j for i, j in enumerate(_nc2code)}
code2ecode = {i : codeToECode(i) for i, j in enumerate(_nc2code)}

# some general helper functions for other program parts
isPair = lambda x, y : code2ecode[x] == y
isPseudoPair = lambda x, y : (code2ecode[x] & y != 0) or (code2ecode[y] & x != 0)
