#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from random import shuffle, choice, randint
from base import _pairs, charsToChar
from Bio.Seq import reverse_complement

"""
Helper stuff for testing.
"""

#
# Random sequence
#
def getNuc(gcContent):
    """
        Get a nuc list.

        @param gcContent The gcContent of the list.
        @result The nuc list.
    """
    nuc = []
    for i in range(100):
        if gcContent * 100 > i: nuc.append('G'); nuc.append('C')
        else:                   nuc.append('A'); nuc.append('T')
    return nuc    

def getRandNucString(gcContent, length):
    """
        Get a random nuc. string.

        @param gcContent The gcContent of the genome/string.
        @param length The length of the nuc. string.
    """
    # make a nuc list corresponding to the gcContent
    nuc = getNuc(gcContent)
    # now construct a random string
    return "".join([choice(nuc) for x in range(length)])    

#
# Sequence evolving
#
def insertSNP(seq, gcContent):
    """
        Insert a single snp into a sequence.

        @param seq The sequence to introduce the snp to.
        @param gcContent The gcContent of the sequence.
        @result The sequence with a SNP.
    """
    i = randint(0, len(seq) - 1)
    return seq[:i] + choice(getNuc(gcContent)) + seq[i+1:]

def insertSNPs(seq, gcContent, nr):
    """
        Let a nuc sequence "evolve" by inserting SNPs.

        @param seq The sequence to evolve.
        @param gcContent The gc content of the genome.
        @param nr The number of SNPs to simmulate.
        @result The evolved sequence.
    """
    for _x_ in range(nr): insertSNP(seq, gcContent)
    return seq

def insertDel(seq, length):
    """
        Let the sequence evolve, by inserting a deletion.

        @param seq The sequence to evolve.
        @param length The length of the region in the sequence to delete.
        @return The evolved sequence.
    """
    if len(seq) <= length: return ''
    i = randint(0, len(seq) - length)
    return seq[:i] + seq[i+length:]

def insertDels(seq, delFreq, nr):
    """
        Let a sequence evolve by inserting deletions into a sequence.

        @param seq The sequence to evolve.
        @param delFreq A list giving the frequence of a deletion length.
        @param nr The number of times to execute the deletion step.
        @result The evolved sequence.
    """
    for _x_ in range(nr): seq = insertDel(seq, choice(delFreq))
    return seq

def insertCnv(seq, length):
    """
        Let the sequence evolve, by inserting a copy.

        @param seq The sequence to evolve.
        @param length The length of the copy.
        @result The evolved sequence.
    """
    i = randint(0, len(seq) - length)
    return seq[:i] + seq[i:i+length] * 2 + seq[i+length:]

def insertCnvs(seq, insFreq, nr):
    """
        Let the sequence evolve, by inserting multiple copies.

        @param seq The sequence to evolve.
        @param insFreq A list giving the frequence of an insertion length.
        @param nr The number of times to execute the insertion step.
        @result The evolved sequence.
    """
    for _x_ in range(nr): seq = insertCnvs(seq, choice(insFreq))
    return seq

#
# misc functions
#
def maybeSwitch(s1, s2):
    """
        Maybe swith sequences.

        @param s1 The first sequence to maybe switch.
        @param s2 The second sequence to maybe switch.
    """
    if choice([0, 1]) == 1:
        return (s1, s2)
    return (s2, s1)

def combineSeqs(s1, s2):
    """
        "Combine" two nuc. sequences.

        @param s1 The first of the two nuc. sequences to combine.
        @param s2 The second of the two nuc. sequences to combine.
        @return The combined sequence.
    """
    r = []
    for i, c in enumerate(s1):
        if c == s2[i]: r.append(c)
        else: r.append(charsToChar([c, s2[i]]))
    return "".join(r)

def actionOnRandom(lst, func):
    """
        Performs an action onto a random element in a list.

        @param string The list of elements.
        @param func The function to perform.
        @result The new list.
    """
    i = randint(0, len(lst) - 1)
    lst[i] = func(lst[i])
    return lst

#
# - calculate a long nuc. sequence (the "reference" genome)
# - get a copy of the reference genome (g1 and g2)
# - let g1 and g2 evolve
# by inserting SNPs and insertions (deletions not nec. because versus visa to insertion for other string)
# when inserting do a starting point correction
# - get sequences for c11, c12 and c21, c22 respectively
# - return data
#
def getTestChromStrings(gcContent, s1, s2, l1, l2, snps, lfreq, cnvs):
    """
        Get two test string representing two read out chromatin strings.

        @param gcContent The gc content of the string.
        @param s1 The starting point for the forward sequence (>0).
        @param s2 The starting point for the reverse sequence (>0).
        @param l1 The length of the forward sequence (>0).
        @param l2 The length of the reverse sequence (>0).
        @param snps The number of snps to insert.
        @param lfreq A list representing the frequence of insertion length.
        @param cnvs The number of cnvs to insert.
        @return The chromatin strings.
    """
    # first step, calculate the reference genome
    refGenome = getRandNucString(gcContent, max(s1+l1, l2)+10+max(lfreq)*cnvs)
    # extract c1, c2 from the reference genome
    g1, g2 = refGenome[:], refGenome[:]
    # evolution
    s11, s12, s21, s22 = s1, s1, s2, s2
    while snps > 0 or cnvs > 0: # as long as the evolution has still some actions
        actions = [] # list with all available actions
        if snps > 0: actions.append(0)
        if cnvs > 0: actions.append(1)
        if choice(actions) == 0:
            g1, g2 = actionOnRandom([g1, g2], lambda x : insertSNP(x, gcContent))
            snps -= 1
        else: # insertion
# insertion don't make sense into the primer regions => should be in "overlapping" regions ...
            if choice([0, 1]) == 0:
                l = choice(lfreq)
                i = randint(max(s11, s21), min(s11+l1, s21+l2) - l)
                g1 = g1[:i] + g1[i:i+l] * 2 + g1[i+l:]
                s21 += l
            else:
                l = choice(lfreq)
                i = randint(max(s12, s22), min(s12+l1, s22+l2) - l)
                g2 = g2[:i] + g2[i:i+l] * 2 + g2[i+l:]
                s22 += l
            cnvs -= 1
    # calc chromatin sequences for forward/reverse + g1/g2
    c11, c12, c21, c22 = g1[s11:s11+l1], g2[s12:s12+l1], g1[s21:s21+l2], g2[s22:s22+l2]
    # combine and return data
    return ((combineSeqs(c11, c12), combineSeqs(c21, c22)), (combineSeqs(reverse_complement(c11), reverse_complement(c12)), combineSeqs(reverse_complement(c21), reverse_complement(c22))), (refGenome, g1, g2), (c11, c12, c21, c22))

if __name__ == "__main__":
    import random
    gcContent = random.random()
    len1 = random.randint(1000, 1500)
    len2 = random.randint(1000, 1500)
    startP1 = random.randint(15, 25)
    startP2 = random.randint(0, 5)
    snps = random.randint(0, 100)
    gap = random.randint(1, 15)
    print(gcContent, len1, len2, startP1, startP2, snps, gap)
    testData = getTestChromStrings(gcContent, startP1, startP2, len1, len2, snps, [gap], 1)
    print(testData)
