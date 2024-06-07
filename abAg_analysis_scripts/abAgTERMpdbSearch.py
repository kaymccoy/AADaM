import sys
import argparse
sys.path.append('/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/lib')
import mstpython as mst



parser = argparse.ArgumentParser()
parser.add_argument("-s", "--struct", dest="struct", help="the path to the structure to analyze - antibody must be IMGT numbered")
parser.add_argument("-abC", "--abChains", dest="abChains", help="the antibody chains, comma separated")
parser.add_argument("-agC", "--agChains", dest="agChains", help="the antigen chains, comma separated")
parser.add_argument("-o", "--outFile", dest="outFile", help="The outFile, to save all info to")
parser.add_argument("-r", "--rmsdRest", dest="r", help="The max RMSD for matches")
arguments = parser.parse_args()



def addOrAppend(v,d):
    if v in d:
        d[v] += 1
    else:
        d[v] = 1



def abAgSearch(setChainsDiffBool,seqRestrictBool,myFASST,abBool):

    nonPolarList = ['G','A','V','C','P','L','I','M','W','F']
    polarList = ['S','T','Y','N','Q']
    plusChargeList = ['K','R','H']
    minusChargeList = ['D','E']

    aaTypeDict = {'G': nonPolarList, 'A': nonPolarList, 'V': nonPolarList,'C' : nonPolarList, 'P': nonPolarList, 'L': nonPolarList, 'I': nonPolarList, 'M': nonPolarList, 'W': nonPolarList, 'F': nonPolarList, 'S': polarList, 'T': polarList, 'Y': polarList, 'N': polarList, 'Q': polarList, 'K': plusChargeList, 'R': plusChargeList, 'H': plusChargeList, 'D': minusChargeList,'E': minusChargeList}

    # load structure / set up chains

    structPath = arguments.struct
    S = mst.Structure(structPath, "QUIET")
    abChains = arguments.abChains.split(',')
    if arguments.agChains:
        agChains = arguments.agChains.split(',')
    else:
        agChains = []
        chainCount = S.chainSize()
        i = 0
        while i < chainCount:
            c = S.getChain(i)
            i += 1
            cName = c.id
            if cName in abChains:
                continue
            else:
                agChains.append(cName)

    # go over each residue on antibody side in crystal (later make this a sub-function and do for model as well as crystal)

    numMatchesList = []
    RMSDsList = []

    abResiduesDict = {}
    agResiduesDict = {}
    IMGTdict = {}
    numContactsList = []
    
    fullAbResiduesDict = {}
    cdrAbResiduesDict = {}
    fullAgResiduesDict = {}
    alreadyCounted = []

    abResIDs = []
    abXYZs = []

    sResList = S.getResidues()
    cf = mst.ConFind("/dartfs/rc/lab/G/Grigoryanlab/home/coy/MST/testfiles/rotlib.bin", S)

    for i in range(len(sResList)):

        numContacts = 0

        iRes = sResList[i]
        try:
            iAtom = iRes.findAtom("CA",True)
        except:
            continue
        abX = iAtom.x
        abY = iAtom.y
        abZ = iAtom.z

        iName = iRes.name
        iChain = iRes.getChainID(True)
        
        if iChain in agChains:
            try:
                addOrAppend(mst.SeqTools.tripleToSingle(iName,"+"),fullAgResiduesDict) 
            except:
                continue

        if iChain not in abChains:
            continue
        
        try:
            addOrAppend(mst.SeqTools.tripleToSingle(iName,"+"),fullAbResiduesDict) 
        except:
            continue

        abNum = iRes.num
        abInsert = iRes.iCode
        if (27 <= abNum and abNum <= 38) or (56 <= abNum and abNum <= 65) or (105 <= abNum and abNum <= 117):
            addOrAppend(mst.SeqTools.tripleToSingle(iName,"+"),cdrAbResiduesDict) 


        cl = mst.ContactList()
        try:
            cf.getResidueContacts(iRes,0.001,cl)
        except:
            continue
        cl.sortByDegree()

        realAbRes = ''
        resList = []
        segMid1 = 1
        abConstList = []
        segMid2 = 1
        agConstList = []
        for ii in range(len(cl)):
            resList = []
            currRes = cl.residueB(ii)
            currChain = currRes.getChainID(True)

            # if in contact w/ antigen chain...

            if currChain in agChains:

                # get the ab side flanking residues

                piFlankingRes = iRes.previousResidue()
                if piFlankingRes != None:
                    resList.append(piFlankingRes)
                else:
                    segMid1 = 0

                resList.append(iRes)
                realAbRes = iRes.name
                abConstList = [mst.SeqTools.tripleToSingle(realAbRes,"+")]

                try:
                    niFlankingRes = iRes.nextResidue()
                    if niFlankingRes != None:
                        resList.append(niFlankingRes)
                except:
                    niFlankingRes = None

                # get the ag side flanking residues

                piiFlankingRes = currRes.previousResidue()
                if piiFlankingRes != None:
                    resList.append(piiFlankingRes)
                else:
                    segMid2 = 0

                resList.append(currRes)
                realAgRes = iRes.name
                agConstList = aaTypeDict[mst.SeqTools.tripleToSingle(realAgRes,"+")]

                try:
                    niiFlankingRes = currRes.nextResidue()
                    if niiFlankingRes != None:
                        resList.append(niiFlankingRes)
                except:
                    niiFlankingRes = None

                break

        # if no matches on ag side found, just continue to the next residue

        if len(resList) == 0:
            continue
        
        for ii in range(len(cl)):
            currRes = cl.residueB(ii)
            currChain = currRes.getChainID(True)
            if currChain in agChains:
                numContacts += 1
                currIndex = currRes.getResidueIndex()
                if currIndex not in alreadyCounted:
                    alreadyCounted.append(currIndex)
                    currAgRes = currRes.name
                    try:
                        addOrAppend(mst.SeqTools.tripleToSingle(currAgRes,"+"),agResiduesDict)
                    except:
                        continue
        numContactsList.append(numContacts)
        addOrAppend(mst.SeqTools.tripleToSingle(realAbRes,"+"),abResiduesDict)

        if (27 <= abNum and abNum <= 38) or (56 <= abNum and abNum <= 65) or (105 <= abNum and abNum <= 117):
            addOrAppend(str(abNum) + str(abInsert),IMGTdict)
        else:
            addOrAppend("non-IMGT",IMGTdict)

        termStruct = mst.Structure(resList)

        try:
            myFASST.setQuery(termStruct,False) #False = don't split by chain according to geometry
        except:
            continue
        
        if setChainsDiffBool:
            myFASST.options.setChainsDiff(0,1)
        else:
            myFASST.options.setMinGap(0,1,1)

        if seqRestrictBool:
            seqConst = mst.fasstSeqConstSimple(myFASST.numQuerySegments)
            seqConst.addConstraint(0, segMid1, abConstList)
            seqConst.addConstraint(1, segMid2, agConstList)
            myFASST.options.setSequenceConstraints(seqConst)

        myFASST.search()

        numMatchesBase = myFASST.numMatches
        abResIDs.append(str(iChain) + ";" + str(abNum) + ";" + str(abInsert))
        abXYZs.append(str(abX) + ";" + str(abY) + ";" + str(abZ))

        myMatches = myFASST.getMatches()
        topMatch = myMatches[0]
        topMatchRMSDbase = topMatch.rmsd

        # NOW SEARCH CLOSEST TERM

        myFASST.options.resetGapConstraints(2)
        myFASST.options.resetDiffChainConstraints(2)

        topMatchStruct = myFASST.getMatchStructure(topMatch, False, mst.matchType.WITHGAPS, False)
        #try:
        myFASST.setQuery(topMatchStruct,False) #False = don't split by chain according to geometry
        #except:
        #    continue

        if setChainsDiffBool:
            myFASST.options.setChainsDiff(0,1)
        else:
            myFASST.options.setMinGap(0,1,1)

        myFASST.search()

        numMatchesSecond = myFASST.numMatches
        numMatchRatio = float(numMatchesSecond+1)/float(numMatchesBase+1)
        numMatchesList.append(numMatchRatio)

        myMatches = myFASST.getMatches()
        topMatch = myMatches[0]
        topMatchRMSDsecond = topMatch.rmsd
        topMatchRatio = float(topMatchRMSDsecond)/float(topMatchRMSDbase)
        RMSDsList.append(topMatchRatio)

    if setChainsDiffBool:
        myFASST.options.resetDiffChainConstraints(2)
    else:
        myFASST.options.resetGapConstraints(2)

    outFileStr = ""

    if abBool:
        outFileStr += "AB-AG "
    else:
        outFileStr += "PDB40 "
    if setChainsDiffBool:
        outFileStr += " CHAINS DIFF RESTRICTED + "
    else:
        outFileStr += " CHAINS NOT RESTRICTED + "
    if seqRestrictBool:
        outFileStr += "SEQ RESTRICTED: "
    else:
        outFileStr += "SEQ NOT RESTRICTED: "
    outFileStr += "number of matches per TERM,"
    for i in numMatchesList:
        outFileStr += str(i) + ","
    outFileStr += "\n"

    if abBool:
        outFileStr += "AB-AG "
    else:
        outFileStr += "PDB40 "
    if setChainsDiffBool:
        outFileStr += " CHAINS DIFF RESTRICTED + "
    else:
        outFileStr += " CHAINS NOT RESTRICTED + "
    if seqRestrictBool:
        outFileStr += "SEQ RESTRICTED: "
    else:
        outFileStr += "SEQ NOT RESTRICTED: "
    outFileStr += "top match RMSD per TERM,"
    for i in RMSDsList:
        outFileStr += str(i) + ","
    outFileStr += "\n"

    extraAnalysis = False
    if setChainsDiffBool: #write extra analysis only once

        outFileStr += "number of contacts between antibody residue and antigen residues,"
        for i in numContactsList:
            outFileStr += str(i) + ","
        outFileStr += "\n"

        outFileStr += "IMGT contact,"
        nextLine = "IMGT contact residue count,"
        for key in IMGTdict:
            outFileStr += str(key) + ","
            nextLine += str(IMGTdict[key]) + ","
        outFileStr += "\n" + nextLine + "\n"

        outFileStr += "antibody contact residues,"
        nextLine = "antibody contact residue count,"
        for key in abResiduesDict:
            outFileStr += str(key) + ","
            nextLine += str(abResiduesDict[key]) + ","
        outFileStr += "\n" + nextLine + "\n"

        outFileStr += "antibody full residues,"
        nextLine = "antibody full residue count,"
        for key in fullAbResiduesDict:
            outFileStr += str(key) + ","
            nextLine += str(fullAbResiduesDict[key]) + ","
        outFileStr += "\n" + nextLine + "\n"

        outFileStr += "CDR full residues,"
        nextLine = "CDR full residue count,"
        for key in cdrAbResiduesDict:
            outFileStr += str(key) + ","
            nextLine += str(cdrAbResiduesDict[key]) + ","
        outFileStr += "\n" + nextLine + "\n"

        outFileStr += "antigen contact residues,"
        nextLine = "antigen contact residue count,"
        for key in agResiduesDict:
            outFileStr += str(key) + ","
            nextLine += str(agResiduesDict[key]) + ","
        outFileStr += "\n" + nextLine + "\n"

        outFileStr += "antigen full residues,"
        nextLine = "antigen full residue count,"
        for key in fullAgResiduesDict:
            outFileStr += str(key) + ","
            nextLine += str(fullAgResiduesDict[key]) + ","
        outFileStr += "\n" + nextLine + "\n"

        outFileStr += "ab res ID info,"
        for i in abResIDs:
            outFileStr += str(i) + ","
        outFileStr += "\n"
        
        outFileStr += "xyz coords,"
        for i in abXYZs:
            outFileStr += str(i) + ","
        outFileStr += "\n"

    return outFileStr



finalOutFile = ''



#searchDB = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/tinyRedundancyTestDb.sim.bin'
searchDB = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/PDB40_fixed.bin'
searchFASST = mst.FASST()
print("reading database...")
searchFASST.readDatabase(searchDB,0)
print("done reading database...")
searchFASST.options.redundancyCut = 0.5
searchFASST.options.rmsdCutoff = float(arguments.r)
searchFASST.options.minNumMatches = 1
searchFASST.options.maxNumMatches = 100



print("searching... chains the same, unrestricted by sequence...")
finalOutFile += abAgSearch(False,False,searchFASST,False)

print("searching... chains different, unrestricted by sequence...")
finalOutFile += abAgSearch(True,False,searchFASST,False)

#print("searching... chains the same, restricted by sequence...")
#finalOutFile += abAgSearch(False,True,searchFASST,False)

#print("searching... chains different, restricted by sequence...")
#finalOutFile += abAgSearch(True,True,searchFASST,False)



with open(arguments.outFile,"w") as of:
        of.write(finalOutFile)
