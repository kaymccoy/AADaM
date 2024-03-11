# code for parsing down antibody-protein structures since a date D, removing those with matches already in SAbDab (and so the PDB) with antigens at X% seq ID or less to their antigen /and/ an antibody at Y% seq ID or less to their antibody

import argparse
import datetime
import glob
import sys
sys.path.append('/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/lib')
import mstpython as mst
from anarci import anarci
from Bio import Align
from Bio.Seq import Seq
import os
import pickle
import shutil


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputDb", dest="inputDb", help="path to the input database, which should be a complete download of all antibodies in complex with protein antigens, and have a summary file ending in 'summary.tsv'")
parser.add_argument("-d", "--date", dest="date0", help="the date on or after which you want to make your test set; format YYYY/MM/DD")
parser.add_argument("-cd", "--cutoffDate", dest="cutoffDate", help="a date past which you won't accept structures")
parser.add_argument("-c1", "--abCompSeqCut", dest="abCompSeqCut", help="percent sequence ID minimum to disallow an ab-ag complex, based on past ab-ag complexes, comparing H loops to H loops and L loops to L loops")
parser.add_argument("-c2", "--withinDatasetCut", dest="withinDatasetCut", help="percent sequence ID minimum used for knocking out complexes that are too similar w/ in the database, based on H or L or any antigen chain to any other antigen chain sequence ID")
parser.add_argument("-m", "--methodsAllowed", dest="methodsAllowed", help="methods allowed, separated by a comma, like 'X-RAY DIFFRACTION,ELECTRON MICROSCOPY'")
parser.add_argument("-w", "--whiteList", dest="whiteList", help="a comma separated list of pdb files that you want to definitely be included in the dataset (will still be filtered by resolution, sequence identity, method, etc. - this just ensures when structures 'knock each other out' the whitelisted ones will be given preference to remain. Useful for adding on new structures to a previously-made dataset)'")
parser.add_argument("-o", "--outDir", dest="outDbPath", help="")
parser.add_argument("-r", "--resCut", dest="resCut", help="defaults to 100")
parser.add_argument("-cs", "--cutoffStrict", dest="cutoffStrict", help="make the cutoffs for sequence ID within the group strict, by applying them to the antigen and antibody loop sequences individually, instead of in combination (i.e. strict cutoffs would not allow an antibody-antigen complex with the same antigen as a previously accepted structure but different antibody in, while standard / not strict cutoffs would allow it)")
parser.add_argument("-nx", "--nx", dest="nx", help="strip unnatural residues from the ends, and don't allow sequences with unnatural residues in the middle")
parser.add_argument("-g", "--globalSeqID", dest="globalSeqID", help="use global sequence ID, number of matches divided by the length of the shorter sequence, to get sequence ID %. This can make the cutoff more stringent in practice, as a shorter sequence can map well to a much larger non-similar one by allowing many gaps. The default (when this flag is not given) is to do local alignment with gap penalites, then apply the sequence identity cutoffs to the best local alignment, if it's over 30 residues long.")
arguments = parser.parse_args()
arguments = parser.parse_args()

date0 = arguments.date0.split("/")
d0 = datetime.datetime(int(date0[0]), int(date0[1]), int(date0[2]))

if arguments.cutoffDate:
    cdList = arguments.cutoffDate.split("/")
    cd = datetime.datetime(int(cdList[0]), int(cdList[1]), int(cdList[2]))

if not os.path.exists(arguments.outDbPath):
    os.makedirs(arguments.outDbPath)

if not arguments.resCut:
    arguments.resCut = 100

if arguments.whiteList:
    whiteList = arguments.whiteList.split(",")
else:
    whiteList = []

percentWinDatasetCut = float(arguments.withinDatasetCut)/100
pAbSeqCut = float(arguments.abCompSeqCut)/100.0

testPrintName = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/abProtParseDownTestWriting.txt'



globalXXaligner = Align.PairwiseAligner()
localMSaligner = Align.PairwiseAligner()
localMSaligner.mode = "local"
localMSaligner.match = 2
localMSaligner.mismatch_score = -1
localMSaligner.open_gap_score = -0.5
localMSaligner.extend_gap_score = -0.1




def seqPer(a,aa,alignType):
    try:
        if alignType == 'local':
            alignmentList = localMSaligner.align(a, aa)
            alignment = alignmentList[0]
            string = alignment._format_psl()
            lst = string.split('\t')
            matchNum = int(lst[0])
            gaps1num = int(lst[5])
            gaps2num = int(lst[7])
            denom = matchNum + gaps1num + gaps2num
            if denom < 30:
                return 0
                
        elif alignType == 'global':
            matchNum = globalXXaligner.score(a, aa)
            denom = min(len(a),len(aa))
        
        else:
            print("can only provide global or local as the alignment type...")
            quit()

        per = float(matchNum/denom)
        return per
    
    except:
        print("no alignment made...")
        print("a was:")
        print(a)
        print("aa was:")
        print(aa)
        return 0.0
        # the above accounts for cases like a: 'RSISSINIHTRDGSTTTLTGFPRIRS', aa: 'AAAAAAAAAAAAAAAAAAAAAAAAAAY' wherein there's no decent local alignment at all... in which case you should just skip checking the sequence ID % and return 0



def isScFv(hChain,lChain):
    hChainCapital = hChain.upper()
    lChainCapital = lChain.upper()
    if hChain == lChainCapital:
        return True
    elif lChain == hChainCapital:
        return True
    else:
        return False



def tripleToSingle(aa):
    if aa == "HIS":
        return "H"
    elif aa == "ALA":
        return "A"
    elif aa == "GLY":
        return "G"
    elif aa == "ILE":
        return "I"
    elif aa == "LEU":
        return "L"
    elif aa == "PRO":
        return "P"
    elif aa == "VAL":
        return "V"
    elif aa == "PHE":
        return "F"
    elif aa == "TRP":
        return "W"
    elif aa == "TYR":
        return "Y"
    elif aa == "ASP":
        return "D"
    elif aa == "GLU":
        return "E"
    elif aa == "ARG":
        return "R"
    elif aa == "LYS":
        return "K"
    elif aa == "SER":
        return "S"
    elif aa == "THR":
        return "T"
    elif aa == "CYS":
        return "C"
    elif aa == "MET":
        return "M"
    elif aa == "ASN":
        return "N"
    elif aa == "GLN":
        return "Q"
    else:
        #print("unknown amino acid " + aa)
        return "X"

# parses sequence from SEQRES section of a PDB file


def seqResParser(filePath):
    seqDict = {}
    pdbFile = open(filePath, 'r')
    fileList = pdbFile.readlines()
    for fileLine in fileList:
        lineSplit = fileLine.split()
        if lineSplit[0] == 'SEQRES':
            lineChain = fileLine[11]
            lineSeq = ''
            lineLen = len(lineSplit)
            a = 4
            while a < lineLen:
                lineRes = lineSplit[a]
                shortRes =tripleToSingle(lineRes)
                lineSeq += shortRes
                a += 1
            if lineChain in seqDict:
                seqDict[lineChain] += lineSeq
            else:
                seqDict[lineChain] = lineSeq
    return seqDict



# anarci parser

def anarciParser(oriSeqRes,anarciOut):

    # go along the tuples of the original sequence, and check the next residue in the anarci output; if it matches, the IMGT numbering to the tuple, and incriment to move on to the next in the original sequence; if it's a -, continue, and if it's neither... report an error cuz something's gone wrong

    # oriI & oriEnd will be values like 0 & 115, where oriI is "The index in the sequence of the first numbered residue", so essentially - when you put in the sequence AMG... if the numbering starts at M, oriI = 1, and oriEnd is "The index in the sequence of the last numbered residue", so if your sequence is 120 residues long but the last thing numbered by anarci is at index 114, oriEnd will be 114.

    # note: as anNumbering (the numbering according to anarci) will sometimes have lots of dashes in it (that is, instead of matching a residue you put in to a number, it has a gap / dash there that matches to an number) it will not necessarily be the size of the span of indexes (like 0-115 = 116), it will be at least that length or higher

    try:
        anNumbering = anarciOut[0]
        oriI = anarciOut[1]
        oriEnd = anarciOut[2]
    except:
        print("oriSeqRes: ")
        print(oriSeqRes)
        print("anNumbering: ")
        print(anNumbering)

    oriSeqIMGTtuples = []

    for anRes in anNumbering:
        anNumber = anRes[0]
        anName = anRes[1]
        if anName == '-':
            continue
        elif anName == oriSeqRes[oriI]: #sanity check to make sure it's the right residue
            tupleToAdd = (anNumber,anName)
            oriSeqIMGTtuples.append(tupleToAdd)
            oriI += 1
        else:
            print("error: anarci numbering & original pdb numbering off somehow")
            print("oriI:")
            print(oriI)
            print("oriEnd")
            print(oriEnd)
            print("anNumbering:")
            print(anNumbering)
            print("anarciOut:")
            print(anarciOut)
            print("oriSeqRes")
            print(oriSeqRes)
            print()
            quit()

    return oriSeqIMGTtuples



# structs to seqs function


def abAgStructs2Seqs(listArg,inputDb):

    infoAndSeqList = []

    for i in listArg:

        currFilePath = inputDb + "/" + str(i[1]) + ".pdb"

        #currS = str(i[0])
        currPDB = str(i[1])
        currAg = str(i[4])
        currAbH = str(i[2])
        currAbL = str(i[3])
        currComplexType = str(i[5])
        currRes = i[6]

        if isScFv(currAbH,currAbL):
            abScFv = True
            currAbH = currAbH.upper()
            currAbL = currAbL.upper()
        else:
            abScFv = False

        # load the struct

        try:
            structAbAg = mst.Structure(currFilePath, "QUIET")
        except:
            print("testing failed file path to get structure from")
            print(currFilePath)
            continue

        seqResDict = seqResParser(currFilePath)

        # 1. get the sequences for the antibody and antigen [antigen seq,HloopSeqs,LloopSeqs] format

        if currAbH == "NA" and currAbL == "NA":
            print("WARNING: skipping entry" + currPDB + " because it apparently has neither heavy or light chains")
            continue

        seqFail = False

        agList = []
        aSeqReses = []
        aSeqs = []
        aPenalty = 0

        if currAg != "NA":
            if " | " in currAg:
                agList = currAg.split(" | ")
            else:
                agList = [currAg]
            for ag in agList:
                aChain = structAbAg.getChainByID(ag)
                aSeq = ''
                try:
                    resA = aChain.getResidues()
                except:
                    seqFail = True
                    print("WARNING: failed at getting antigen residues from entry " + currPDB + "with chain ID(s): ")
                    print(agList)
                    continue

                try:
                    aSeqRes = seqResDict[ag]
                except:
                    seqFail = True
                    print("WARNING: failed at getting antigen sequence from entry " + currPDB + "with ag chain: " + ag)
                    print("dictionary of all sequences is as follows: ")
                    print(seqResDict)
                    continue

                #aNames = []
                for ii in range(len(resA)):

                    aRes = resA[ii]
                    aName = aRes.name
                    aNameS = tripleToSingle(aName)
                    aSeq += aNameS

                aSeqReses.append(aSeqRes)
                aSeqs.append(aSeq)

            a = 0
            aPenalty = 0
            while a < len(aSeqReses):
                aPenaltyScore = seqPer(aSeqs[a], aSeqReses[a],'global')          
                minAlen = min(len(aSeqs[a]), len(aSeqReses[a]))
                aPenalty += aPenaltyScore - minAlen
                a +=1

        else:
            print("WARNING: no antigen found for entry " + currPDB)

        if seqFail == True:
            continue

        if currAbH != "NA":
            chainAbH = structAbAg.getChainByID(currAbH)
            try:
                resH = chainAbH.getResidues()
            except:
                print("WARNING: getting heavy chain residues failed for entry " + currPDB + " with chain ID: " + currAbH)
                continue
            hNames = []
            hSeq = ''
            if currAbH not in seqResDict:
                #print("chains off in file; skipping")
                continue
            hSeqRes = seqResDict[currAbH]

            for ii in range(len(resH)):
                hRes = resH[ii]
                hName = hRes.name
                hNum = hRes.num
                hCode = hRes.iCode
                hIndex = hRes.getResidueIndex()
                hNameS =tripleToSingle(hName)
                hNames.append([hIndex,(hNum,hCode),hNameS])
                hSeq += hNameS

            hPenaltyScore = seqPer(hSeq,hSeqRes,'global')
            minHlen = min(len(hSeq), len(hSeqRes))
            hPenalty = hPenaltyScore - minHlen
        else:
            hSeqRes = ''
            hPenalty = 0

        if abScFv:
            lSeqRes = hSeqRes
            lPenalty = 0
        elif currAbL != "NA":
            chainAbL = structAbAg.getChainByID(currAbL)
            try:
                resL = chainAbL.getResidues()
            except:
                print("WARNING: getting light chain residues failed for entry " + currPDB + " with chain ID: " + currAbL)
                continue
            lNames = []
            lSeq = ''

            if currAbL not in seqResDict:
                continue

            lSeqRes = seqResDict[currAbL]

            for ii in range(len(resL)):
                lRes = resL[ii]
                lName = lRes.name
                lNum = lRes.num
                lCode = lRes.iCode
                lIndex = lRes.getResidueIndex()
                lNameS =tripleToSingle(lName)
                lNames.append([lIndex,(lNum,lCode),lNameS])
                lSeq += lNameS
            lPenaltyScore = seqPer(lSeq, lSeqRes,'global')
            minLlen = min(len(lSeq), len(lSeqRes))
            lPenalty = lPenaltyScore - minLlen
        else:
            lSeqRes = ''
            lPenalty = 0

        # 2. run anarci to get alignment, and 3. correspond the old numbering & indexes to Anarci numbering so loops can be identified later

        hNumbering = "NA"
        if currAbH != "NA":
            anarciH, anarciHextra, _ = anarci([("placeholderName",hSeqRes)], scheme="imgt", output=False)

            anarciHdicts = anarciHextra[0]
            if anarciHdicts is None:
                #print('anarciHdicts empty; skipping')
                continue

            d = 0
            anarciHl = 'NA'
            while d < len(anarciHdicts):
                anarciDict = anarciHdicts[d]
                if "chain_type" in anarciDict and anarciDict['chain_type'] == 'H':
                    anarciHl = anarciH[0][d]
                    break
                else:
                    d += 1

            if anarciHl == 'NA':
                print("H chain in " + currPDB +"'s anarci output not found, this structure will be skipped'")
                print("anarciHdicts were:")
                print(anarciHdicts)
                print("from h chain sequence:")
                print(hSeqRes)
                continue

            hNumbering = anarciParser(hSeqRes,anarciHl) # removed sorted() around anarciParser output

        lNumbering = "NA"
        if currAbL != "NA":
            anarciL, anarciLextra, _ = anarci([("placeholderName",lSeqRes)], scheme="imgt", output=False)
            anarciLdicts = anarciLextra[0]

            if anarciLdicts is None:
                continue

            d = 0
            anarciLl = 'NA'
            while d < len(anarciLdicts):
                anarciDict = anarciLdicts[d]
                if "chain_type" in anarciDict and ((anarciDict['chain_type'] == 'L') or (anarciDict['chain_type'] == 'K')):
                    anarciLl = anarciL[0][d]
                    break
                else:
                    d += 1

            if anarciLl == 'NA':
                print("L chain in " + currPDB +"'s anarci output not found, this structure will be skipped'")
                print("anarciLdicts were:")
                print(anarciLdicts)
                print("from l chain sequence:")
                print(lSeqRes)
                continue

            lNumbering = anarciParser(lSeqRes,anarciLl) # removed sorted() around anarciParser output

        # parse down the sequences to only be IMGT numbered ones

        hSeqShort = ''
        lSeqShort = ''
        if hNumbering != 'NA':
            for h in hNumbering:
                if len(h) == 1:
                    continue
                else:
                    hSeqShort += h[1]

        if lNumbering != 'NA':
            for l in lNumbering:
                if len(l) == 1:
                    continue
                else:
                    lSeqShort += l[1]

        # get all the loops by IMGT; CDR1 = 27-38, CDR2 = 56-65, CDR3 = 105-117

        currCDRH1 = ''
        currCDRH2 = ''
        currCDRH3 = ''
        currCDRL1 = ''
        currCDRL2 = ''
        currCDRL3 = ''
        if currAbH != "NA":
            ni = 27
            while ni < 118:
                for hTuple in hNumbering:
                    if len(hTuple) == 2:
                        hIMGTnum = hTuple[0][0]
                        if float(hIMGTnum) == ni:
                            if 27 <= float(hIMGTnum) and float(hIMGTnum) <= 38:
                                currCDRH1 += hTuple[1]
                            if 56 <= float(hIMGTnum) and float(hIMGTnum) <= 65:
                                currCDRH2 += hTuple[1]
                            if 105 <= float(hIMGTnum) and float(hIMGTnum) <= 117:
                                currCDRH3 += hTuple[1]
                ni += 1
            hLoopSeqs = currCDRH1 + currCDRH2 + currCDRH3
            hLoopSeqs = hLoopSeqs.strip("X")
        else:
            hLoopSeqs = ''

        if currAbL != "NA":
            ni = 27
            while ni < 118:
                for lTuple in lNumbering:
                    if len(lTuple) == 2:
                        lIMGTnum = ''
                        lIMGTnum = lTuple[0][0]
                        if float(lIMGTnum) == i:
                            if 27 <= float(lIMGTnum) and float(lIMGTnum) <= 38:
                                currCDRL1 += lTuple[1]
                            if 56 <= float(lIMGTnum) and float(lIMGTnum) <= 65:
                                currCDRL2 += lTuple[1]
                            if 105 <= float(lIMGTnum) and float(lIMGTnum) <= 117:
                                currCDRL3 += lTuple[1]
                ni += 1
            lLoopSeqs = currCDRL1 + currCDRL2 + currCDRL3
            lLoopSeqs = lLoopSeqs.strip("X")
        else:
            lLoopSeqs = ''

        totScore = aPenalty + hPenalty + lPenalty

        if arguments.nx:
            newaSeqs = []
            for aSeq in aSeqReses:
                newSeq = aSeq.strip("X")
                if newSeq != '':
                    newaSeqs.append(newSeq)
            hSeqShort = hSeqShort.strip("X")
            lSeqShort = lSeqShort.strip("X")

        infoAndSeqList.append({'currPDB':currPDB, 'aSeq':newaSeqs, 'hSeq': hSeqShort, 'lSeq': lSeqShort, 'hLoopSeqs':hLoopSeqs, 'lLoopSeqs':lLoopSeqs, 'complexType':currComplexType, 'currAg':currAg, 'currAbH':currAbH, 'currAbL':currAbL, 'totScore':totScore, 'currRes':currRes})

    return infoAndSeqList

cutoffFilterCount = 0
methodFilterCount = 0
qualFilterCount = 0
prevSimilarityCount = 0
knockoutWithinCount = 0

# load in summary file

onOrAfterDateList = []
finalList = []
beforeDateProtPepComplexesList = []

methodsList = arguments.methodsAllowed.split(",")

abAgSysIndex = -1
sumFilePath = arguments.inputDb + "/*.tsv"
print(sumFilePath)
sumFileName = glob.glob(sumFilePath)
if len(sumFileName) !=1:
    print("\n" + "You have the wrong number of summary files, there should only be 1. Check your abag dir")
with open(sumFileName[0]) as file:
    unsplitLines = file.read().splitlines()
    for unsplitLine in unsplitLines:
        abAgSysIndex += 1
        line = unsplitLine.split("\t")
        currPdb = line[0]
        if currPdb == 'pdb':
            continue
        if len(line) < 18:
            continue
        
        currAbH = line[1];
        currAbL = line[2];
        currAg = line[4];
        currComplexType = line[5];
        currDate = line[9];

        dateSplit = currDate.split('/')
        dateYearEnd = int(dateSplit[2])
        if dateYearEnd > 50:
            dateYear = dateYearEnd + 1900
        else:
            dateYear = dateYearEnd + 2000
        try:
            d1 = datetime.datetime(dateYear, int(dateSplit[0]), int(dateSplit[1]))
        except:
            #print("date out of range so skipping line:")
            #print(line)
            continue
        
        try:
            currResolution = float(line[16]);
        except:
            qualFilterCount += 1
            continue
        
        currMethod = line[17];

        # make a list A of all the PDB IDs on or after D

        currList = [currDate,currPdb,currAbH,currAbL,currAg,currComplexType,currResolution]

        if arguments.cutoffDate:
            if d1 >= cd:
                continue

        if d1 >= d0:
            
            if currMethod not in methodsList:
                methodFilterCount += 1
                continue

            if currResolution > float(arguments.resCut):
                #print("skipping " + currPdb + " because resolution too poor")
                qualFilterCount += 1
                continue
            if 'protein' in currComplexType:
                onOrAfterDateList.append(currList)

        # make a list B of all the PDB IDs before D

        else:
            if 'protein' in currComplexType:
                cutoffFilterCount += 1
            if 'protein' in currComplexType or 'peptide' in currComplexType:
                beforeDateProtPepComplexesList.append(currList)

print("cutoffFilterCount:")
print(cutoffFilterCount)
print("qualFilterCount:")
print(qualFilterCount)
print("methodFilterCount:")
print(methodFilterCount)

# get sequences etc. from pdb files - each entry gets a dict

onOrAfterDateInfoAndSeqsList = abAgStructs2Seqs(onOrAfterDateList,arguments.inputDb)

print("onOrAfterDateInfoAndSeqsList len")
print(len(onOrAfterDateInfoAndSeqsList))

if len(beforeDateProtPepComplexesList) > 0:
    beforeDatePpInfoAndSeqsList = abAgStructs2Seqs(beforeDateProtPepComplexesList,arguments.inputDb)
else:
    beforeDatePpInfoAndSeqsList = []

print("beforeDatePpInfoAndSeqsList len")
print(len(beforeDatePpInfoAndSeqsList))

acceptedAbProtStructs = []

for i in onOrAfterDateInfoAndSeqsList:
    toRemoveIfAdding = []
    keepI = True
    acceptedStructInfo = {}

    if arguments.nx:
        if ("X" in i['hSeq']) or ("X" in i['lSeq']) or ("X" in i['aSeq']):
            continue

    currPdb = i['currPDB']

    print("working on classifying entry: " + currPdb)

    ChSeq = Seq(i['hLoopSeqs'])
    ClSeq = Seq(i['lLoopSeqs'])
    CaSeqs = [Seq(s) for s in i['aSeq']]
    currNegScore = i['totScore']
    currComplexType = i['complexType']
    currRes = i['currRes']

    # check for repeats
    worseRepeat = False
    betterRepeat = False
    ii = 0

    while ii < len(acceptedAbProtStructs):

        # when considering which structures to add, compare H loop sequences, L loop seqeunces, and antigen sequences to previously accepted structures

        closestAbMatchScore = 0
        closestAgMatchScore = 0

        Ppdb = acceptedAbProtStructs[ii]['currPDB']
        PhSeq = Seq(acceptedAbProtStructs[ii]['hLoopSeqs'])
        PlSeq = Seq(acceptedAbProtStructs[ii]['lLoopSeqs'])
        PaSeqs = [Seq(s) for s in acceptedAbProtStructs[ii]['aSeq']]
        pastNegScore = acceptedAbProtStructs[ii]['totScore']
        pastRes = acceptedAbProtStructs[ii]['currRes']

        if len(PhSeq) > 1 and len(ChSeq) > 1:
            if arguments.globalSeqID:
                hScore = seqPer(PhSeq,ChSeq,'global')
            else:
                hScore = seqPer(PhSeq,ChSeq,'local')

            if closestAbMatchScore < hScore:
                closestAbMatchScore = hScore

        if len(PlSeq) > 1 and len(ClSeq) > 1:
            if arguments.globalSeqID:
                lScore = seqPer(PlSeq,ClSeq,'global')
            else:
                lScore = seqPer(PlSeq,ClSeq,'local')

            if closestAbMatchScore < lScore:
                closestAbMatchScore = lScore


        aSimBest = 0
        caLen = 0


        for a in CaSeqs:
            if len(a) == 0:
                continue
            caLen += len(a)
            paLen = 0
            for aa in PaSeqs:
                paLen += len(aa)

                if arguments.globalSeqID:
                    aScore = seqPer(a,aa,'global')
                else:
                    aScore = seqPer(a,aa,'local')
                
                if closestAgMatchScore < aScore:
                    closestAgMatchScore = aScore
        
        if (arguments.cutoffStrict and ((closestAgMatchScore >= percentWinDatasetCut) or (closestAbMatchScore >= percentWinDatasetCut))) or (not arguments.cutoffStrict and ((closestAgMatchScore >= percentWinDatasetCut) and (closestAbMatchScore >= percentWinDatasetCut))):

            knockoutWithinCount += 1

            # if something previous that it might knockout is in the whitelist, don't add it and keep the whitelisted structure. Then if the prev structure isn't whitelisted, keep whichever one has the fewest breaks in the H / L loops and antigen; if equal amounts of breaks then keep the one with the shortest antigen; if equal antigen lengths, keep the previously accepted one!
            
            if Ppdb in whiteList:
                worseRepeat = True
                break

            if currNegScore < pastNegScore:
                worseRepeat = True
                break
            elif currNegScore > pastNegScore:
                toRemoveIfAdding.append(ii)
                if len(toRemoveIfAdding) > 1:
                    worseRepeat = True
                    break
                else:
                    betterRepeat = True
                    ii += 1
                    continue # if a better repeat may not be the best so keep checking
            else:
                if paLen <= caLen:
                    worseRepeat = True
                    break
                else:
                    toRemoveIfAdding.append(ii)
                    if len(toRemoveIfAdding) > 1:
                        worseRepeat = True
                        break
                    else:
                        betterRepeat = True
                        ii += 1
                        continue 
        else:
            ii += 1

    if worseRepeat == True:
        continue

    for ii in beforeDatePpInfoAndSeqsList:
        PhSeq = Seq(ii['hLoopSeqs'])
        PlSeq = Seq(ii['lLoopSeqs'])
        prevPDB = ii['currPDB']
        hAlighments = []
        bestHalighn = 0.0
        lAlighments = []
        bestLalighn = 0.0

        # check if the sequence for the antibody is at X% seq ID or less to a past example in the before list - so has this antibody's "binding mode" been seen before?

        pAbScore = 0.0

        if len(i['hLoopSeqs']) > 0 and len(ii['hLoopSeqs']) > 0:

            if arguments.globalSeqID:
                pAbScore = seqPer(ChSeq,PhSeq,'global')
            else:
                pAbScore = seqPer(ChSeq,PhSeq,'local')

        if len(i['lLoopSeqs']) > 0 and len(ii['lLoopSeqs']) > 0:

            if arguments.globalSeqID: 
                bestLalign = seqPer(ClSeq,PlSeq,'global')
            else:
                bestLalign = seqPer(ClSeq,PlSeq,'local')

            if bestLalign > pAbScore:
                pAbScore = bestLalign

        if (pAbScore >= pAbSeqCut):
            prevSimilarityCount
            keepI = False
            break

    # if everything's checked and there're no matches, start an accepted structure entry

    if keepI == True:

        acceptedStructInfo = i

        if len(toRemoveIfAdding) > 1: # if you'd be knocking out multiple structures to only include one more, that's NG
            continue

        acceptedAbProtStructs.append(acceptedStructInfo)

        # remove past worst duplicate, if relevant

        if len(toRemoveIfAdding) > 0:
            acceptedAbProtStructsNew = []
            iii = 0
            while iii < len(acceptedAbProtStructs):
                if iii not in toRemoveIfAdding:
                    acceptedAbProtStructsNew.append(acceptedAbProtStructs[iii])
                iii += 1
            acceptedAbProtStructs = acceptedAbProtStructsNew

print("prevSimilarityCount: ")
print(prevSimilarityCount)
print("knockoutWithinCount: ")
print(knockoutWithinCount)

# save & print the accepted structures

print("done! writing files...")

with open(arguments.outDbPath + "/fullDb.pkl", 'wb') as f:
    pickle.dump(acceptedAbProtStructs, f)

f = open(arguments.outDbPath + "/lightDb.txt", 'w')
f.write("pdb ID,A chain(s),H chain(s),L chain(s)\n")
for i in acceptedAbProtStructs:
    f.write(i['currPDB'] + "," + i['currAg'] + "," + i['currAbH'] + "," + i['currAbL'] + "\n")

structPath = arguments.outDbPath + "/structures"
if not os.path.exists(structPath):
    os.makedirs(structPath)
for i in acceptedAbProtStructs:
    currFilePath = arguments.inputDb + "/" + str(i['currPDB']) + ".pdb"
    shutil.copyfile(currFilePath, structPath + "/" + str(i['currPDB']) + ".pdb")

# write complex fastas

if not os.path.exists(arguments.outDbPath + "/complexFastas"):
    os.makedirs(arguments.outDbPath + "/complexFastas")

for i in acceptedAbProtStructs:
    fileName = arguments.outDbPath + "/complexFastas/" + i['currPDB'] + "complex.fasta"
    f = open(fileName, 'w')
    count = 1
    for ag in i['aSeq']:
        f.write("> Antigen " + str(count) + "\n")
        f.write(ag + "\n")
        count += 1
    if i['hSeq'] == i['lSeq']:
        f.write("> Antibody H and L seq\n")
        f.write(i['hSeq'] + "\n")
    else:
        if i['hSeq'] != '':
            f.write("> Antibody H\n")
            f.write(i['hSeq'] + "\n")
        if i['lSeq'] != '':
            f.write("> Antibody L\n")
            f.write(i['lSeq'] + "\n")
    f.close()

# write complex fastas colon seperated

compColSepPath = arguments.outDbPath + "/complexFastasColonSep"
if not os.path.exists(compColSepPath):
    os.makedirs(compColSepPath)

for i in acceptedAbProtStructs:
    fileName = compColSepPath + "/" + i['currPDB'] + "complexColSep.fasta"
    f = open(fileName, 'w')
    f.write("> Complex \n")
    startBool = True
    for ag in i['aSeq']:
        if startBool:
            f.write(ag)
            startBool = False
        else:
            f.write(":" + ag)

    if i['hSeq'] == i['lSeq']:
        f.write(":" + i['hSeq'])
    else:
        if i['hSeq'] != '':
            f.write(":" + i['hSeq'])
        if i['lSeq'] != '':
            f.write(":" + i['lSeq'])
    f.close()

# write antigen fastas

if not os.path.exists(arguments.outDbPath + "/antigenFastas"):
    os.makedirs(arguments.outDbPath + "/antigenFastas")

for i in acceptedAbProtStructs:
    fileName = arguments.outDbPath + "/antigenFastas/" + i['currPDB'] + "antigen.fasta"
    f = open(fileName, 'w')
    count = 1
    for ag in i['aSeq']:
        f.write("> Antigen " + str(count) + "\n")
        f.write(ag + "\n")
        count += 1
    f.close()

# write antigen fastas colon seperated

antiColSepPath = arguments.outDbPath + "/antigenFastasColonSep"
if not os.path.exists(antiColSepPath):
    os.makedirs(antiColSepPath)

for i in acceptedAbProtStructs:
    fileName = antiColSepPath + "/" + i['currPDB'] + "antigenColSep.fasta"
    f = open(fileName, 'w')
    f.write("> Antigen \n")
    startBool = True
    for ag in i['aSeq']:
        if startBool:
            f.write(ag)
            startBool = False
        else:
            f.write(":" + ag)
    f.close()

# write antibody fastas

if not os.path.exists(arguments.outDbPath + "/antibodyFastas"):
    os.makedirs(arguments.outDbPath + "/antibodyFastas")

for i in acceptedAbProtStructs:
    fileName = arguments.outDbPath + "/antibodyFastas/" + i['currPDB'] + "antibody.fasta"
    f = open(fileName, 'w')
    if i['hSeq'] == i['lSeq']:
        f.write("> Antibody H and L\n")
        f.write(i['hSeq'] + "\n")
    else:
        if i['hSeq'] != '':
            f.write("> Antibody H\n")
            f.write(i['hSeq'] + "\n")
        if i['lSeq'] != '':
            f.write("> Antibody L\n")
            f.write(i['lSeq'] + "\n")
    f.close()

# write antibody fastas colon seperated

abColSepPath = arguments.outDbPath + "/antibodyFastasColonSep"
if not os.path.exists(abColSepPath):
    os.makedirs(abColSepPath)

for i in acceptedAbProtStructs:
    fileName = abColSepPath + "/" + i['currPDB'] + "antibodyColSep.fasta"
    f = open(fileName, 'w')
    f.write("> Antibody \n")
    startBool = True
    if i['hSeq'] == i['lSeq']:
        f.write(i['hSeq'])
    else:
        if i['hSeq'] != '':
            f.write(i['hSeq'])
            startBool = False
        if i['lSeq'] != '':
            if startBool:
                f.write(i['lSeq'])
            else:
                f.write(":" + i['lSeq'])
    f.close()

# print number of structures

print("base number of structure past date: " + str(len(onOrAfterDateInfoAndSeqsList)))
print("total number of structures accepted at end: " + str(len(acceptedAbProtStructs)))
