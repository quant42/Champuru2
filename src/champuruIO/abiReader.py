#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
from Bio._py3k import _bytes_to_string, _as_bytes
from Bio import SeqIO

__doc__ = """
Reader implementation for the abi file format specification.
"""

def isReadable(handle):
    """
        Fast check, if this file is readable by this reader. Check if the
        file magic bytes equals to 'ABIF' as specified in the file format
        specification:
        http://www6.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf

        @param handle The file handle.
        @return True if this (probably) is an abi file.
    """
    # check if first 4 bytes equals to ABIF
    return handle.read(4) == _as_bytes('ABIF')

def readFile(handle):
    """
        Read in the abi file.

        @param handle The file handle.
        @return The file data.
    """
    # use Biopython in order to read in the abif file.
    record = SeqIO.read(handle, "abi")
    # convert to interface
    # get "Sequencing Analysis Filter wheel order"
    # should be GATC, but check anyway
    nns = record.annotations['abif_raw']['FWO_1']
    # transform
    data = {
        "FILE_EXT" : "ABI",
        "DATA" : record
    }
    for i, char in enumerate(nns):
        if char not in ['A', 'C', 'T', 'G']:
            raise Exception("Unexpected character in file '%s'!" % char)
        data[char] = record.annotations['abif_raw']['DATA' + str(9 + i)]
    # return the data
    return data

