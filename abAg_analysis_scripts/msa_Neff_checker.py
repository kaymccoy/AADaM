# looks at average Neff for single chain MSA, total Neff for single chain MSA, and Neff for paired-MSA

import glob
import matplotlib.pyplot as plt
#from Bio import pairwise2
#import string
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdbID", dest="pdbID", help="pdbID to run on")
arguments = parser.parse_args()               

# tricky side note: insertions change the "line-up" of each sequence to the other, so they need to be removed - but counted, as a mismatch. That's why the lower case characters (insertions) are removed and added to the currCount

def getNeff(list):

    newDic = {}
    prevKeys = []
    seqCount = 101
    seqKey = 'NA'
    li = 0
    while li < len(list):

        # if the expected next / new key...
        
        if (list[li].strip() == ">" + str(seqCount)) or (list[li].strip() == "\x00>" + str(seqCount)):
            #print("testing 1...")
            seqKey = str(seqCount) # set seqKey here
            prevKeys.append(seqKey)
            seqCount += 1 # so makes 101 become 102 etc.
            li += 2 # skip this entry b/c it's the original sequence rather than a match
            continue
        
        # if a repeat of a previous key...

        elif (list[li][0] == ">") or ("\x00" in list[li][0]):
            for p in prevKeys:
               foundPrev = False
               if (list[li].strip() == ">" + str(p)) or (list[li].strip() == "\x00>" + str(p)):
                #print("testing 2...")
                seqKey = p # set seqKey here
                seqCount += 1
                li += 2 # skip its seq cuz not a match, it's the original sequence
                foundPrev = True
                continue
            if not foundPrev: # if it's a match (because it's not a new or old key)
                li += 1
                continue
        
        else: # now we're on the match sequence
            if seqKey in newDic:
                newDic[seqKey].append(list[li].strip())
            else:
                newDic[seqKey] = [list[li].strip()]

        li += 1

    # go over every line and calc nef for each

    neffList = []

    for key in newDic:

        neffCurr = 0

        l = newDic[key]
        print("len for key: " + key)
        print(len(l))
        i = 1
        while i < len(l):

            if i % 1000 == 1:
                print("another 1000 lines done...")

            simCount = 0

            currSeq = l[i]

            # go over every other line to compare to

            ii = 1
            while ii < len(l):

                compSeq = l[ii].strip()

                currPer = 0
                currCount = 0
                currTot = 0
                iii = 0
                iv = 0
                
                while (iii < len(currSeq)):

                    currR = currSeq[iii]
                    compR = compSeq[iv]

                    if (currR == "-") and (compR == '-'):
                        iii += 1
                        iv += 1
                        continue
                    
                    elif (currR.islower()) and (not compR.islower()):
                        currTot += 1
                        iii += 1
                        continue
                    
                    elif (currR.islower()) and (not compR.islower()):
                        currTot += 1
                        iv += 1
                        continue   

                    elif currR == compR:
                        currCount += 1
                    
                    currTot += 1
                    iii += 1
                    iv += 1

                if currTot > 0:
                    currPer = float(currCount) / float(currTot)
                else:
                    currPer = 0
            
                # if sequence similarity is >= 0.80, add to total - PREVIOUSLY DID 0.62!
                
                if currPer >= 0.80:
                    simCount += 1

                ii += 1 # iterate over comparision
            
            if simCount > 0: 
                neffCurr += 1 / simCount

            i += 1 # iterate over line to check

        neffList.append(neffCurr)

    if len(neffList) > 0: 
        neffAvg = sum(neffList)/len(neffList)
    else:
        neffAvg = 0.0
    return(neffAvg)

rdpPath = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/AlphaFoldAnalysis/' + arguments.pdbID + '_1_*_RLRDP.txt'
dqPath = '/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/AlphaFoldAnalysis/' + arguments.pdbID + '_1_RDP.txt'
rdpFile = glob.glob(rdpPath)[0]
dqFile = glob.glob(dqPath)[0]

pdbID = arguments.pdbID

# get DOCKQ score / RDP-val

with open(dqFile,"r") as rf:
    dq = float(rf.readlines()[0].strip())

with open(rdpFile,"r") as rf:
    rdp = float(rf.readlines()[1].strip())

# go through all MSAs + get numbers of entries for each (# lines minus 1, then divided by 2) - POSSIBLY READ AB/AG ONLY STUFF LATER AS WELL

scUnirefMSA = glob.glob('/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/AlphaFoldRawOut/complexOut*/' + pdbID + 'complexColSep_env/uniref.a3m')[0]
with open(scUnirefMSA,"r") as su:
    sLines = su.readlines()

scBfdMSA = glob.glob('/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/AlphaFoldRawOut/complexOut*/' + pdbID + 'complexColSep_env/bfd.mgnify30.metaeuk30.smag30.a3m')[0]
with open(scBfdMSA,"r") as sb:
    sbLines = sb.readlines()   
sLines.extend(sbLines)

#print(scBfdMSA)

pcMSA = glob.glob('/dartfs/rc/lab/G/Grigoryanlab/home/coy/AbAgComplexResults/AlphaFoldRawOut/complexOut*/' + pdbID + 'complexColSep_/pair.a3m')[0]
with open(pcMSA,"r") as p:
    pLines = p.readlines()   

#print(pcMSA)

# get n_eff for each

sNeff = getNeff(sLines)
pNeff = getNeff(pLines)

print(dq)
print(rdp)
print(sNeff)
print(pNeff)