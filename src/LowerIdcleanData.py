import re
import os

pathOfSequence="..//NewData//fasta(1).txt"

if(0):
    patt = r"^.....|PDBID$"
    text = list()
    ProteinSeq = list()
    ProteinName = list()
    for line in open(pathOfSequence):
        text.append(line)

    sequence = ""
    name = ""
    num = 0
    for line in text:
        if line.startswith('>'):
            ProteinSeq.append(sequence)
            ProteinName.append(line)
            sequence=''
            num += 1
        else:
            sequence += line

    ProteinSeq.remove(ProteinSeq[0])
    ProteinName.remove(ProteinName[-1])
    # add all names and sequences into lists


    print(num)
    print(len(ProteinSeq))
    print(len(ProteinName))
    for i in range(len(ProteinName)):
        ProteinName[i] = re.findall(patt,ProteinName[i])[0]
        ProteinName[i] = ProteinName[i].strip().replace('>','')
        ProteinName[i] = ProteinName[i].replace('|PDBID|CHAIN|SEQUENCE', '')
        ProteinSeq[i] = ProteinSeq[i].strip().replace('\n', '')
    print(ProteinName)

    i=0
    while i<len(ProteinName):
        if ProteinName[i] == ProteinName[i-1]:
            ProteinName.remove(ProteinName[i])
            ProteinSeq.remove(ProteinSeq[i])
            i -= 1
        i += 1
    print(ProteinName)
    print(len(ProteinSeq))
    print(len(ProteinName))

    path = "..//NewData//NewSeq.data"
    f = open(path,'a')
    for i in range(len(ProteinName)):
        line = ProteinName[i]+','+ProteinSeq[i]+'\n'
        f.write(line)
    f.close()

if(0):
    pair = r'[0-9][A-Z][]'
    Disulpath='..//NewData//SSBond//proteins/'
    files=os.listdir(Disulpath)
    fileinfo=list()
    for file in files:
        bondinfo = list()
        realPath=Disulpath+file
        filename=file[0:-4]
        for line in open(realPath,'r'):
            if line.startswith('SSBOND'):
                line=line[0:50]
                splitPattern=r'CYS .'
                #item=re.split(splitPattern,line)
                #assert len(item)==3  test the status of dataset
                bondinfo.append(line)
        fileinfo.append((filename,bondinfo))





    print(fileinfo)


