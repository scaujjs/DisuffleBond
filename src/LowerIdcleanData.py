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
    files=os.listdir('../PDB2')
    for file in files:
        print(file)
        if ('.ent') in file:
            DSSPname=file[3:7]
        else:
            DSSPname=file[0:-4]
        print(DSSPname)
        prody.execDSSP("../PDB2/"+file,'../DSSP2/'+DSSPname)


listOfpidWithChainId=list()
for line in open('../listOfCulled'):
    listOfpidWithChainId.append(line.strip())

patternSS='NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)'
patternStart='RESIDUE AA STRUCTURE'
def pareDSSPfile(file):
    examples=list()
    start=0
    seq=''
    ###################
    pid=file[0:-5]
    ################
    chainId=''
    path="../DSSP2/"+file
    ##print(path)

    for line in open(path):
        if patternStart in line:
            start=1
            continue
        if start:
            aa=line[13:15].strip()

            if aa!='!*':
                assert len(aa)==1
            if aa=='!*':
                examples.append((pid+':'+chainId,seq))
                seq=''
            else:
                chainId=line[10:13].strip().lower()
                seq+=aa
    examples.append((pid+':'+chainId,seq))
    newExamples=list()
    ##print('helloworld')
    ##print(examples)
    for pidcid,seq in examples:
        if pidcid in listOfpidWithChainId:
            newExamples.append((pidcid,seq))
            print(pidcid,seq)
    return newExamples


if 0:
    files=os.listdir("../DSSP2/")
    lines=list()
    num=0
    for file in files:
        ##for line in open("../DSSP/"+file,'r'):
        exampleTemp=pareDSSPfile(file)
        assert len(exampleTemp)!=0
        for pidcid,seq in exampleTemp:
            strTemp=pidcid+','+seq+'\n'
            lines.append(strTemp)
            num+=1

        ##seqs.append(pareDSSPfile(file))
    print(num)

    f=open('../NewData/seqFromDSSPculled25','w')
    f.writelines(lines)
    f.close()

if 1:
    numOfE=0
    p=r'[a-z]'
    records=list()
    bonds=list()
    newSeq=list()
    keys=list()

    for line in open('../NewData/seqFromDSSPculled25','r'):
        key,seq=line.strip().split(',')
        dictPair=dict()
        odd=0
        for aa in seq:
            if aa=='!':
                continue
            if re.match(p,aa):
                if aa in dictPair:
                    dictPair[aa]=dictPair[aa]+1
                else:
                    dictPair[aa]=1

        for lowcase,value in dictPair.items():
            if value!=2:
                odd=1
        s=seq.replace("!",'')
        if not odd:
            bond=list()
            for lowcase, value in dictPair.items():
                pair=list()
                for j in range(len(s)):
                    if s[j]==lowcase:
                        pair.append(j)
                assert len(pair)==2
                bond.append((pair[0]+1,pair[1]+1))
            newS=s
            for lowcase, value in dictPair.items():
                newS=newS.replace(lowcase,'C')
            bonds.append(bond)
            newSeq.append(newS)
            keys.append(key)

        assert len(newSeq)==len(bonds)


    for i in range(len(newSeq)):
        records.append((keys[i],newSeq[i],bonds[i]))

    print(len(records))

    if 1:
        import random
        random.shuffle(records)
        newSeqs=list()
        newBonds=list()
        for i in range(len(records)):
            newSeqsStr=records[i][0]+','+records[i][1]+'\n'
            bondStr=records[i][0]+','
            print(bondStr)
            for p1,p2 in records[i][2]:
                pairStr=str(p1)+'_'+str(p2)
                bondStr+=(pairStr+',')
            bondStr=bondStr[0:-1]+'\n'
            newSeqs.append(newSeqsStr)
            newBonds.append(bondStr)

        f=open("../NewData/newSeqFromPDBculled25",'w')
        f.writelines(newSeqs)
        f.close()
        f=open("../NewData/newBondsFromPDBculled25",'w')
        f.writelines(newBonds)
        f.close()

def loadTrainDataDSSP():
    listOfKey=list()
    listOfSeq=list()
    listOfPairs=list()

    for line in open("../NewData/newSeqFromDSSP"):
        key,seq=line.strip().split(',')
        listOfSeq.append(seq)
        listOfKey.append(key)
    for line in open("../NewData/newBondsFromDSSP"):
        items=line.strip().split(",")
        pairs=list()
        assert len(items)>0
        for i in range(len(items)-1):
            p1,p2=items[i+1].split("_")
            pairs.append((int(p1),int(p2)))
        listOfPairs.append(pairs)

    return(listOfKey,listOfSeq,listOfPairs)

listOfKey,listOfSeq,listOfPairs=loadTrainDataDSSP()

if 0:
    typeOfStructure=list()
    newlistOfKey = list()
    newlistOfSeq = list()
    newlistOfPairs = list()
    for i in range(len(listOfSeq)):
        for p1,p2 in listOfPairs[i]:
            assert listOfSeq[i][p1-1]=='C'
            assert listOfSeq[i][p2-1]=='C'
        if listOfKey[i][0:4] not in typeOfStructure:
            typeOfStructure.append(listOfKey[i][0:4])
            newlistOfKey.append(listOfKey[i])
            newlistOfSeq.append(listOfSeq[i])
            newlistOfPairs.append(listOfPairs[i])
    print(len(newlistOfPairs))
    recordOfSeqlist=list()
    recordOfBonds=list()
    for i in range(len(newlistOfPairs)):
        recordStrOfSeq=newlistOfKey[i]+','+newlistOfSeq[i]+'\n'
        recordStrOfBond=newlistOfKey[i]+','
        for p1,p2 in newlistOfPairs[i]:
            recordStrOfBond+=(str(p1)+"_"+str(p2)+',')
        recordStrOfBond=recordStrOfBond[0:-1]+'\n'
        recordOfSeqlist.append(recordStrOfSeq)
        recordOfBonds.append(recordStrOfBond)
    assert len(recordOfSeqlist)==len(recordOfBonds)


    f=open("../NewData/newSeqFromPDB_unique",'w')
    f.writelines(recordOfSeqlist)
    f.close()
    f=open("../NewData/newBondsFromPDB_unique",'w')
    f.writelines(recordOfBonds)
    f.close()



    print(len(typeOfStructure))


if 0:
    import tarfile
    import shutil

    filesOfRa = os.listdir('../ra/')
    filesOfPSB2 = os.listdir('../PDB2/')
    for file in filesOfPSB2:
        if '.gz' in file:
            PDBid=file[0:4]
            if PDBid+'.pdb' in filesOfPSB2 or 'pdb'+PDBid+'.ent' in filesOfPSB2:
                pass
            else:
                print(file)

                command='gzip -d ../PDB2/'+file
                os.system(command)

