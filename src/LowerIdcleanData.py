import re
import os
pathOfSequence="../NewData/fasta(1).txt"

if(0):
    patt = r">.*|PDBID"
    text = list()
    ProteinSeq = list()
    ProteinName = list()
    for line in open(pathOfSequence):
        text.append(line)

    sequence = ""
    name = ""
    num = 0
    isFrist=1
    for line in text:
        if re.match(patt,line):
            print(line)
            assert ':' in line
            assert '|PDBID|CHAIN|SEQUENCE' in line
            newLine=line.replace('>','')
            newLine=newLine.replace('|PDBID|CHAIN|SEQUENCE','')
            pid,key=newLine.split(':')
            pid=pid.strip()
            key=key.strip()
            ProteinName.append((pid,key))
            if isFrist:
                isFrist=0
            else:
                ProteinSeq.append(sequence)
                sequence=''
            num += 1
        else:
            templine = line.strip()
            sequence+=templine
    ProteinSeq.append(sequence)
    print(len(ProteinName),len(ProteinSeq))


    path = "..//NewData//NewSeq.data"
    f = open(path,'w')
    for i in range(len(ProteinName)):
        id=ProteinName[i][0]+":"+ProteinName[i][1]
        line =id+','+ProteinSeq[i]+'\n'
        f.write(line)
    f.close()

if(0):
    Disulpath='..//PDB/'
    files=os.listdir(Disulpath)
    fileinfo=list()
    for file in files:
        assert len(file)==8
        errorFlag=0
        bondinfo = list()
        realPath=Disulpath+file
        filename=file[0:-4]
        ssbondInfo=list()
        findSSbond=0
        for line in open(realPath,'r'):
            if line.startswith('SSBOND'):
                findSSbond=1
                newLine=line.strip().split(" ")
                items=list()
                for item in newLine:
                    if item.strip()=='':
                        pass
                    else:
                        items.append(item.strip())


                assert re.match(r"[0-9]*",items[1])
                assert (re.match(r"[0-9]*",items[4]) and re.match(r"[0-9]*",items[7]))
                ssbondInfo.append(items)
            else:
                if findSSbond:
                    break
        bondChainT=list()
        for items in ssbondInfo:
            if items[3] not in bondChainT:
                bondChainT.append(items[3])
            if items[6] not in bondChainT:
                bondChainT.append(items[6])
        trueComp=list()
        for t in bondChainT:
            errorFlag=0
            for items in ssbondInfo:
                if items[3]==t and items[6]!=t:
                    errorFlag=1
                if items[3]!=t and items[6]==t:
                    errorFlag=1
            if not errorFlag:
                trueComp.append(t)
        infos=list()
        for t in trueComp:
            infoString=''
            for items in ssbondInfo:
                if items[3]==t:
                    info=items[4]+'_'+items[7]+','
                    infoString+=info
            infos.append(infoString[:-1])


        assert len(infos)==len(trueComp)

        onlyTrueType=list()
        onlyInfo=list()
        for i in range(len(infos)):
            if infos[i] not in onlyInfo:
                onlyTrueType.append(trueComp[i])
                onlyInfo.append(infos[i])

        onlyTrueTypeOnlyCYS=list()
        onlyInfoOnlyCYS=list()
        for t in range(len(onlyTrueType)):
            errorFlag=0
            for items in ssbondInfo:
                if items[3]==onlyTrueType[t]:
                    if items[2]!='CYS':
                        errorFlag=1
                    if items[5]!='CYS':
                        errorFlag=1
            if not errorFlag:
                onlyTrueTypeOnlyCYS.append(onlyTrueType[t])
                onlyInfoOnlyCYS.append(onlyInfo[t])





        for i in range(len(onlyTrueTypeOnlyCYS)):
            newRecord=file[0:4]+':'+onlyTrueTypeOnlyCYS[i]+","
            newRecord+=onlyInfoOnlyCYS[i].strip()+'\n'
            fileinfo.append(newRecord)
            print(fileinfo[-1])
    print(len(fileinfo))

    f=open('../NewData/newBond','w')
    f.writelines(fileinfo)
    f.close()

if 0:
    ##clear self connect, clear a cys in more bond
    newBondv1=list()
    for line in open('../NewData/newBond','r'):
        items=line.split(",")
        positionShowBefore=list()
        errorFlag=0
        for i in range(len(items)-1):
            p1,p2=items[i+1].split("_")
            if p1 not in positionShowBefore:
                positionShowBefore.append(p1)
            else:
                errorFlag=1
            if p2 not in positionShowBefore:
                positionShowBefore.append(p2)
            else:
                errorFlag=1

        if not errorFlag:
            newBondv1.append(line)
        else:
            print(line)

    print(len(newBondv1))
    f=open("../NewData/newBondv1",'w')
    f.writelines(newBondv1)
    f.close()





##detect specialPosition error
if 0:
    newLines=list()
    for line in open("../NewData/newBondv1",'r'):
        items=line.strip().split(",")
        file,seqChaim=items[0].split(":")
        pairs=list()
        errorFlag=0
        for i in range(len(items)-1):
            p1,p2=items[i+1].split("_")
            try:
                p1=int(p1)
                p2=int(p2)
            except ValueError:
                errorFlag=1
            ##pairs.append((p1,p2))
        if not errorFlag:
            newLines.append(line)
        else:
            print(line)
    f=open('../NewData/newBondv2','w')
    f.writelines(newLines)
    f.close()
        ##listOfBond.append((file,seqChaim,pairs))
        ##print(listOfBond[-1])
if 0:
    bondDict=dict()
    seqDict = dict()

    newLines=list()
    for line in open("../NewData/newBondv2",'r'):
        items=line.strip().split(",")
        file,seqChaim=items[0].split(":")
        pairs=list()
        for i in range(len(items)-1):
            p1,p2=items[i+1].split("_")
            p1=int(p1)
            p2=int(p2)
            pairs.append((p1,p2))
        key=file+":"+seqChaim
        bondDict[key]=pairs




    for line in open("../NewData/NewSeq.data"):
        items=line.strip().split(',')
        seq=items[1]
        file,seqChaim=items[0].split(":")
        file=file.lower()
        key=file+":"+seqChaim
        seqDict[key]=seq

    keys=list()
    pairs=list()
    seq=list()
    newKey=list()
    for key, pair in bondDict.items():
        keys.append(key)
        pairs.append(pair)
        if key not in seqDict:
            print('error')
        seq.append(seqDict[key])
    print(len(keys),len(pairs),len(seq))
    print(seq[-1])

    num=0
    for i in range(len(keys)):
        errorFlag=0
        seqs=seq[i]
        pairsS=pairs[i]
        for p1,p2 in pairsS:
            p1=int(p1)
            p2=int(p2)
            if 0<=p1<=len(seqs) and 0<=p2<=len(seqs):
                if seqs[p1-1]=='C' and seqs[p2-1]=='C':
                    pass
                else:
                    errorFlag = 1
                    break
            else:
                errorFlag=1
                break
        if errorFlag:
            num+=1
            print(keys[i])
        else:
            newKey.append(keys[i])

    print(len(newKey))
    print(num)




'''
            if seqS[p1-1]!='C':
                print(seqS)
                print(p1)
                print(key)
            if seqS[p2-1]!='C':
                print(seqS)
                print(p2)
                print(key)

'''






if 0:
    numOfLongerThan600=0
    for line in open("../NewData/finalSeq"):
        if len(line)>600:
            numOfLongerThan600+=1
    print(numOfLongerThan600)


AAtrype=list()
def parse2getSeqAndBond(file):
    Disulpath='..//PDB/'
    keys=list()
    seqs=list()
    bonds=list()
    findSeq = 0
    seqFinished=0
    seqTemp = list()
    for line in open(Disulpath+file,'r'):
        if line.startswith('SEQRES   1'):
            if findSeq:
                ##print(len(seqs))

                seqs.append(seqTemp)
                ##print(len(seqs))
            else:
                findSeq=1
            keyC=line[10:12].strip()
            assert keyC not in keys
            assert len(keyC)==1
            keys.append(keyC)
            seqStr=line[19:].strip().split(" ")
            seqTemp=seqTemp+seqStr
        elif line.startswith('SEQRES'):
            seqStr=line[19:].strip().split(" ")
            seqTemp=seqTemp+seqStr
        else:
            if findSeq and (1-seqFinished):
                seqs.append(seqTemp)
                seqFinished=1


    if 0:
        print(file)
        for key in keys:
            print("\t"+key)

    if 1:
        for i in range(len(seqs)):
            for AA in seqs[i]:
                if AA not in AAtrype:
                    print(AA,file,keys[i])
                    AAtrype.append(AA)



    ##print(len(keys))
    ###print(len(seqs))



    return keys


if 0:
    Disulpath='..//PDB/'
    files=os.listdir(Disulpath)
    fileinfo=list()

    numOfTotal=0
    for file in files:
        listOfItems=parse2getSeqAndBond(file)
        numOfTotal+=len(listOfItems)
    print(numOfTotal)
if 0:
    import prody
    files=os.listdir('../PDB')
    for file in files:
        print(file)
        prody.execDSSP("../PDB/"+file,'../DSSP/'+file[0:-4])

patternSS='NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)'
patternStart='RESIDUE AA STRUCTURE'
def pareDSSPfile(path):
    found=0
    numOfChains=0
    numOfTotalSSbond=0
    numOfIntraBond=0
    numOfInterBond=0
    numOfTotalAA=0
    tempChianId='A'
    tempChain=''
    start=0
    for line in open(path):
        ##print(line)
        if patternSS in line:
            items=line.strip().split(" ")
            itemUseful=list()
            for item in items:
                itemTemp=item.strip()
                if itemTemp!='':
                    itemUseful.append()
            numOfChains=int(itemUseful[1])
            numOfTotalSSbond=int(itemUseful[2])
            numOfIntraBond=int(itemUseful[3])
            numOfInterBond=int(itemUseful[4])
            numOfTotalAA=int(itemUseful[0])
        if patternStart in line:
            start=1
        if start:
            itemUseful=list()
            items=line.split(" ")
            chainID=''
            index=''
            AA=''
            for item in items:
                itemTemp=item.strip()
                if itemTemp!='':
                    itemUseful.append(itemTemp)


            chainID=itemUseful[2]
            index=itemUseful[1]
            AA=itemUseful[3]




    return found

files=os.listdir("../DSSP/")
keys=list()
seqs=list()
bonds=list()
for file in files:
    ##print(file)

    key,seq,bond=pareDSSPfile("../DSSP/"+file)
    keys.append(key)
    seqs.append(seq)
    bonds.append(bond)

