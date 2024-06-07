# imports

import argparse
import glob
import sys
import os
import time
sys.path.append('/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/lib') #put the path to your Mosaist lib here!
import mstpython as mst
from anarci import anarci
import string



# arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="i", help='the input file(s) to clean up; if mutiple provide a pattern like "/dartfs-hpc/rc/home/4/f002v94/labhome/AbAgComplexResults/ClusProCFabCFagOut/*/model.000.00_cleaned.pdb", for example')
parser.add_argument("-o", dest="o", help='the output directory to save the outFiles')
parser.add_argument("-n", dest="n", type=str, choices=["dir","file"], default='file', help="the naming pattern for the outFile: either name after the directory they're from (use value 'dir') or the file they're from (use value 'file')")
parser.add_argument("-r", dest="r", type=bool, default=False, help='rename chains (H &/or L for antibody, other chain IDs for the antigen); defaults to False')
parser.add_argument("-s", dest="s", type=bool, default=False, help='snugDock mode - name all antigen chains as A actually, and name in A_HL order, as well as noting if H or L or both chains in final name; should use with -r flag unlelss the H&L chains are already named H/L')
parser.add_argument("-ps", dest="pdbDict", type=string, default=False, help="a set of pdb IDs and their chains to skip renaming (e.g. the antigen is also an antibody but you don't want to rename it); format is pdbID1,chain1,chain2;pdbID2,chain1,chain2;...")
parser.add_argument("-sum", dest="summaryList", default=False, help='if provided, search through a summary list, and only do this process to those files named in the summary file')
arguments = parser.parse_args()

if arguments.pdbDict:
    #split the pdb string by ; and then by , to get the pdbDict
    pdbDict = {}
    for p in arguments.pdbDict.split(";"):
        pSplit = p.split(",")
        pdbDict[pSplit[0]] = pSplit[1:]

# check if out-dir exists; if not, make it

if not os.path.exists(arguments.o):
    os.mkdir(arguments.o)

# chain letters

alphabetSoup = string.ascii_uppercase
alphabetSoupList = list(alphabetSoup)

if arguments.summaryList:
    with open(arguments.summaryList,"r") as sf:
        sumList = sf.readlines()
        pdbList = [x.split(",")[0] for x in sumList]

# read in each structure

fileList = glob.glob(arguments.i)

for f in fileList:
    pdbID = f.split("/")[-1].split("_")[0]

    if arguments.summaryList:
        if pdbID not in pdbList:
            continue

    agCount = 0
    hAlreadyFound = False 
    lAlreadyFound = False

    pdbNameBase = ''
    if arguments.n == 'dir':
        pdbNameBase = f.split("/")[-2] + '_' + f.split("/")[-1]
    else:
        pdbNameBase = f.split("/")[-1]

    print("loading structure for: " + f)

    sAbAg = mst.Structure(f, "QUIET")

    # go through all chains to get all chain IDs:

    c = 0
    chDict = {}
    while c < sAbAg.chainSize():
        ch = sAbAg.getChain(c)
        chID = ch.id
        chDict[chID] = ch
        c += 1

    # go back through all chains: 

    c = 0
    while c < sAbAg.chainSize():
        ch = sAbAg.getChain(c)
        chID = ch.id

        print('working on chain with ID: ' + chID)

        chRs = ch.getResidues()
        chSeq = ''
        chTupList = []
        for i in range(len(chRs)):
            res = chRs[i]
            rName = res.name
            try:
                rSname = mst.SeqTools.tripleToSingle(rName,"+")
            except:
                rSname = 'X'
            if rSname == '?':
                rSname = 'X'
            chSeq += rSname

        # run through anarci to see if it's an antibody

        anarciOut, anarciOutExtra, _ = anarci([("placeholderName",chSeq)], scheme="imgt", output=False)

        if (anarciOut[0] != None):

            if if pdbID in pdbDict and chID in pdbDict[pdbID]:
                c+=1
                continue

            # label the antibody chains if asked to

            if arguments.r:
                if anarciOutExtra[0][0]['chain_type'] == 'H':

                    if hAlreadyFound:
                        print("h chain already found???")
                        quit()
                    hAlreadyFound = True

                    if ('H' in chDict) and (chID != "H"): #if another chain is named H, gotta rename it...
                        pastHchain = chDict["H"]
                        pastRenamed = False
                        for a in alphabetSoup:
                            currLetter = a
                            if a not in chDict: # need to rename to something unique / not previously used
                                pastHchain.id = a
                                chDict[a] = chDict["H"]
                                del chDict["H"]
                                pastRenamed = True
                                break
                        if pastRenamed == False:
                            print("another chain previously named H; unable to rename as there are no other unique letterings available")

                    ch.id = 'H'

                else:

                    if lAlreadyFound:
                        print("l chain already found???")
                        quit()
                    lAlreadyFound = True

                    if ('L' in chDict) and (chID != "L"): #if another chain is named H, gotta rename it...
                        pastLchain = chDict["L"]
                        pastRenamed = False
                        for a in alphabetSoup:
                            currLetter = a
                            if a not in chDict: # need to rename to something unique / not previously used
                                pastLchain.id = a
                                chDict[a] = chDict["L"]
                                del chDict["L"]
                                pastRenamed = True
                                break
                        if pastRenamed == False:
                            print("another chain previously named L; unable to rename as there are no other unique letterings available")

                    ch.id = 'L'

            # parse the anarci output further

            anarciBase = anarciOut[0][0]
            anarciNumbering = anarciBase[0]
            anStart = anarciBase[1]
            anEnd = anarciBase[2]

            # set previous residue numbering to negative

            prevRes = chRs[anStart].previousResidue()
            negCount = -1
            while prevRes != None:
                prevRes.num = negCount
                negCount -= 1
                prevRes = prevRes.previousResidue()

            # set the remainder to IMGT numbering

            for anRes in anarciNumbering:
                anNumber = anRes[0][0]
                anIcode = anRes[0][1]
                anName = anRes[1]

                structRes = chRs[anStart]
                strucResLongName = structRes.name
                srShortName = mst.SeqTools.tripleToSingle(strucResLongName,"+")
                if srShortName == "?":
                    srShortName = 'X'

                if anName == '-':
                    continue
                elif anName == srShortName: #sanity check to make sure it's the right residue
                    structRes.num = int(anNumber)
                    structRes.iCode = anIcode
                    anStart += 1
                else:
                    print("error: anarci numbering & original pdb numbering off somehow")
                    print(anarciOut)
                    print("chSeq: " + chSeq)
                    quit()

            # set further residue numbering

            nextRes = chRs[anEnd].nextResidue()
            nexResNum = chRs[anEnd].num + 1
            while nextRes != None:
                nextRes.num = nexResNum
                nexResNum += 1
                nextRes = nextRes.nextResidue()
        c += 1

    if arguments.r and not hAlreadyFound: #this script intended for nanobodies and heavy + light chain antibodies only; if you want light-chain-only structures as well, mod code as necessary!
        print("no h chain found???")
        quit()

    # fix to snugDock standards if necessary - put all of the antigen chains into one chain A, and have in the LH_A order required by snugdock

    if arguments.s:

        outPathBase = arguments.o + "/" + pdbNameBase
        nS = mst.emptyStructure()
        cc = 0
        rCount = 0
        hChain = None
        lChain = None
        aChain = mst.emptyChain()
        aChain.id = "A"

        while cc < sAbAg.chainSize():
            chch = sAbAg.getChain(cc)
            chchID = chch.id

            if (chchID != "H") and (chchID != "L"):
                chchReses = chch.getResidues()
                r = 0
                while r < len(chchReses):
                    currRes = chchReses[r]
                    currRes.num = rCount
                    rCount += 1
                    r += 1
                    aChain.appendResidue(currRes)

            elif (chchID == "H"):
                hChain = chch
            elif (chchID == "L"):
                lChain = chch
            cc += 1

        if not hChain:
            nS.appendChain(lChain)
            nS.appendChain(aChain)
            outPathBase += "L_A"
        elif not lChain:
            nS.appendChain(hChain)
            nS.appendChain(aChain)
            outPathBase += "H_A"
        else:
            nS.appendChain(lChain)
            nS.appendChain(hChain)
            nS.appendChain(aChain)
            outPathBase += "LH_A"

        outPath = outPathBase + "fixed.pdb"
        nS.writePDB(outPath,"")
    else:
        # save re-numbered & re-chained structure in the out dir
        outPath = arguments.o + "/" + pdbNameBase + "fixed.pdb"
        sAbAg.writePDB(outPath,"")