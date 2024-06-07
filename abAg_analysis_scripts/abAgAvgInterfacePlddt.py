import sys
import argparse
sys.path.append('/dartfs/rc/lab/G/Grigoryanlab/home/coy/Mosaist/lib')
import mstpython as mst



parser = argparse.ArgumentParser()
parser.add_argument("-s", "--struct", dest="struct", help="the path to the structure to analyze - antibody must be IMGT numbered")
parser.add_argument("-abC", "--abChains", dest="abChains", help="the antibody chains, comma separated")
parser.add_argument("-agC", "--agChains", dest="agChains", help="the antigen chains, comma separated")
parser.add_argument("-o", "--outFile", dest="outFile", help="The outFile, to save all info to")
arguments = parser.parse_args()



def addOrAppend(v,d):
    if v in d:
        d[v] += 1
    else:
        d[v] = 1



fullAgResiduesDict = {}
fullAbResiduesDict = {}

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

abAgPlddtL = []
agCounted = []

sResList = S.getResidues()
cf = mst.ConFind("/dartfs/rc/lab/G/Grigoryanlab/home/coy/MST/testfiles/rotlib.bin", S)

for i in range(len(sResList)):

    numContacts = 0

    iRes = sResList[i]
    iChain = iRes.getChainID(True)
    if iChain not in abChains:
        continue
    
    try:
        iAtom = iRes.findAtom("CA",True)
    except:
        continue

    cl = mst.ContactList()
    try:
        cf.getResidueContacts(iRes,0.001,cl)
    except:
        continue
    cl.sortByDegree()

    abCounted = False
    for ii in range(len(cl)):
        currRes = cl.residueB(ii)
        try:
            currAtom = currRes.findAtom("CA",True)
        except:
            continue
        currChain = currRes.getChainID(True)

        # if in contact w/ antigen chain...

        if currChain in agChains:

            if not abCounted:
                abAgPlddtL.append(iAtom.B)
                print("testing iAtom.B: " + str(iAtom.B))

            if currRes in agCounted:
                continue
            else:
                agCounted.append(currRes)
                abAgPlddtL.append(currAtom.B)

print(abAgPlddtL)
averagePLDDT = sum(abAgPlddtL)/float(len(abAgPlddtL))
with open(arguments.outFile,"w") as of:
        of.write(str(averagePLDDT))

# example command: python3 abAgAvgInterPlddt.py -s /dartfs/rc/lab/G/Grigoryanlab/home/coy/databases/abAg_SAbDab_AADaM/structuresIMGT/7amp.pdb -abC H,L -agC A,B -o /dartfs/rc/lab/G/Grigoryanlab/home/coy/Dartmouth_PhD_Repo/plddtInterTest.txt
