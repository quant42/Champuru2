#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import base
from random import gauss
from simSeq import simulate as simSeqSimulate

#seq = "GGCGGCGATCGAGCCGGMWCGCGTCGCCCGCGTGAAGGGTGCGCCGCGCKTCGGGGCGATCACCCGCCGACCGCGCCGCGCGGCKGCGCTTGCCTCTCGCAGCGCCGCTGGGCTGCCCGTGAGGCGCGCSCSGGGTCGCGGCGCTCCTGCGGATGGGCAGGSCGGCTGGAGAGCKGCGGGGTGCGCTCCGCGGGAGACAGGGGCAGTGCGACGGCGCTCGGCCCWGGGGACCGGCCGCTTGGGACCCACKAGGCCCTCAGCGCCTAASCMTCTASGCTYGCTCCGGCGTCTGGGGGCGCGACGAGCGCCCGTCAGGCGGGGGCGCCAAGGTGCGGTAGGGGCTTGGGGCRSCCGTCGTMWTKCGAGGCCGSAGGTGCCCCARACCCCTGGATCGCGCAGTGCCAGCGCTCTTCCGKATCGGCACCGTCCCCCCGACGGCCTACGCCGCGSGCTCGCCCCGGGCGGGCCGCCGCGGCCGCCCCGGGGCGRCACGCCCCCCCCCCCCGGGGGCGCCGSGGATCCCCCCSCGGCCCGCACCGCCGCCGGTGCGGMSGRSSKSSWCGACRSSSTCSWSACSSKSYKGRCGSRGYKSKCGGGGSTSWCKGRGSKSRCYSRRSSSGSSSKRSCGWSGGKSMCGKCSSSCMSGSCCSCYMRGGCCSCKRRKSYSSSKKRTCWSKGCKRKSRMKGCRSKRSMSGCAYARSSCGCAKRKGSCSCRCGYKCCGYSSGCYCSGKCGSSSSYSCSSGSGSSSSSCRSGCSKSCSRGKSYKCSSSGKSYKYSWCGGCMYYCKYSGSRCCSSWSGSSSSGCWSSSSGSRCYCMSGSTWCSCRCSCYTMGSRYCSYKWGSRYSRKSSGMRCGWGCCSRSSGWRSSCGSCKWMGKSGCCYKCKKSGCMSRSGSCSSRSRSGCSCGRSRSKCGCSRGMSYSSSSRGSMCGMGGWGGRCGMGSTGGSYGCSCYSRCWGCYYCCRCWKCYKCCGSSGWSSMMYMRKMGGCGSSSSGSSSGCSSRGYGGYGGSKGGKSGRSSMKCMKMMSRCSSYGCKKSSKGCGKCGYSGSSSGMCSWCSYCSYCSCMRCRRCGRSGRSKKMKKMGYCSMMSRMGSMGSMSGCSSCSCSSCRGRRSSKCCGSCSSSCCSSCCKCYKCYKSCSRCRRYMSSMSYRGKCRGYMSSMMRRARGMSGYCSYYCSKMSGMCSCYCCKMCSRSCSGSSGGGSKGSWSSMSCCSSCGSCKSYYGKCKGSWGGWRTYMKMSSMSSSS"
#seq = "KGGKSGGCGRYSRWSSMGSMRSSYSKCGYCSSCSKSRWGRRKGSKSCGCSSSKCKKSGSGRYSAYCMSCCGMCSRCSSCGCSSSGCKGCGSYKSYYKCYYSYMGCRSCGCYGSKSKGCYSSYSRKGMGSSSSSCSSGKSGYSGCGSYSCTSCKGMKGRKSRGSRSGSCKGSWGRRSWGCGGSGKGSKSYSCKCSGSRGRMRRSRGSRGYRSKRCGRCGSYSSKCSSMSSKGRSSRSCSGCYKSKKRSSMMCYASGMSSYCMKCRSCKMMKMMYCWWCKMKKSYYSCKSCGKCKKSKGGSGSSRCGASSRSCSSYCRKSMGGSGGSGSCRMSRWGSKGYRGKRGSKKSKKGSRGCSSYCGTCKTWMKRSGMSGCMGGWGSYSCMSMMMMCYSSWKSRYSSMGYRSYRSCRSYSYTCYKSMKYRKCRSCRYCSYCCCSMCGRCSKMCKMCGCSGSSYSSYCSCSSSSGGSSSGCCGCSGCSGCCSCSSSGSGGCRMSMCSCCCCCCCCCSSSGGSGSCGSSGMKSMYCCCCCSGCSSSCMSCRCCGCCGSYGSKGMSGSGSWCGRSSWCGRSMYCGWCMYSRYSCKGWSGSSGKRGSYSKSGSKSKYGRSKGWGMSYGMSCGMGCSSGRGYGMSGKYGSSSWCGGSCMCSCSSSCCSSKRMSGCYMMGGSYCYGGGTYWKSRCTRYSRSKGMGGSRSWRSCMSMAYMSGSMKMRGGYGCMCSKSSCKCSSKCKSCSYSGGCSCKGSCSSGSGCSCGSGCSMGGSCYMSGSSTSYGSGYSKGCSWYSKCSKYMYCSGYMSSSSYSSRGSCSSRGSYCSGGSYMCSSRYRCYYMSSSYYCSYASSSYYRKGKKCRYGRGCMSRRSSYGMSCGYGSMMGSCGRMKKCGYSKYCYKSGSCSSMGRSGCMSRGSSCSSGCGMGMGYSRSMGKCSCSSGSMSGMSGAGGWSGRSGYGSKSGCKSYSSCYRCMSCYWCMGCTTCCGCCGKGGACCCATAGGCGGCGCGCGGGCCGMGCGGTGGCGGGTGGGCGAGCCTCAGACCGCGCTGCGTGCGGCGTCGCGGGCCGACCYCGTCCCCGCAACGGCGAGGGCTTAGKCGCCCAAGGCGCAGGCCGCGCCCCGGMAGGKCCGCCGGCCCGCCCGCCTCTGCCYGCGACAGKCCCAGTSGGYAGCCCGAAAGAGGCCGTCCCTCGGACGCCCCTCCGACCGGCCGGGGGGCTGGACCCGCCSGCGCCTGTCGGCTGGAGKTATCCGAGCCCGGGGCTCCTGTCCGGTCCCCTTGGGYGRYGCSGCCG"

def sim(seq):
    peaks = {}
    for pos, char in enumerate(seq):
        chrs = base.codeToChars(base.charToCode(char))
        for chr in chrs:
            peakPos = int(pos * 12 + 5 + gauss(0, 5 / 3) + 0.5)
            lst = peaks.get(chr, [])
            lst.append(peakPos)
            peaks[chr] = lst
    return peaks

if __name__ == "__main__":
    gcContent, len1, len2, startP1, startP2, snps, gap, testdata = simSeqSimulate()
    seqs = testdata[0]
    print(sim(seqs[0]))
    print(sim(seqs[1]))
    print(gcContent, len1, len2, startP1, startP2, snps, gap, testdata)
