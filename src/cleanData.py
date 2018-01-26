import re




pathOfProV="..//data//protVec_100d_3grams.csv"
pathOfBond="..//data//bond.tab"
pathOfSequence="..//data//sequene.fasta"
p1 = r"(?<=<h1>).+?(?=<h1>)"




##this block clean protein sequence
if(0):
    p = r"(?<=\|).+?(?=\|)"
    text=list()
    proteinSeq=list()
    proteinName=list()
    for line in open(pathOfSequence):
        text.append(line)

    sequence=''
    name=''
    num=0
    for line in text:
        if line.startswith('>sp'):
            proteinSeq.append(sequence)
            proteinName.append(line)
            sequence=''
            num = num + 1
        else:
            sequence+=line

    proteinSeq.remove(proteinSeq[0])
    proteinName.remove(proteinName[-1])


    print num
    print len(proteinSeq)
    print len(proteinName)
    for i in range(len(proteinName)):
        proteinName[i]=re.findall(p,proteinName[i])[0]
        ##print proteinName[i]
        proteinSeq[i]=proteinSeq[i].strip().replace('\n','')
        ##print proteinSeq[i]
        ##print

    path="..//data//sequene.data"
    f=open(path,'a')
    for i in range(len(proteinName)):
        line=proteinName[i]+','+proteinSeq[i]+'\n'
        f.write(line)
    f.close()

'''

'''




## this block clean Disulfide bond
if(0):
    pair = r"DISULFID [0-9]+ [0-9]+"
    dis=r"DISULFID"
    lines=list()
    for line in open(pathOfBond):
        lines.append(line)
    lines.remove(lines[0])


    if(1):
        num = 0
        f = open("..//data//bond.data", 'a')
        for line in lines:
            if len(re.findall(dis, line)) == len(re.findall(pair, line))&len(re.findall(dis, line))!=0:
                line2write = ''
                name = line.split()[0]
                line2write = name
                bonds = re.findall(pair, line)
                for b in bonds:
                    line2write = line2write + ',' + b.replace('DISULFID', '')
                line2write = line2write + '\n'
                num = num + 1
                print line2write
                f.write(line2write)
        f.close()
        print num
"""
    for i in range(20):
        print len(re.findall(pair,lines[i]))
        print len(re.findall(dis,lines[i]))
        print lines[i]
"""

##this block is to used clean sequence
if(0):

    index=0
    name=list()
    num=0
    for line in open("..//data//bond.data"):
        name.append(line.split(',')[0])
##        print line.split(',')[0]
        num+=1
    ##print num
    newFasta=list()
    for line in open('..//data//sequene.data'):
        identify=line.split(',')[0]
        if identify==name[index]:
            newFasta.append(line)
            index+=1

    f= open('..//data//cleanSequene.data','a')
    for line in newFasta:
        f.write(line)

    f.close()


##this block is to statistic the num of bond
##all 32674
##max 166
##moreThan99 32
##moreThan9 1613
##moreThan4 6560
##moreThan3 10122
##moreTHan2 16087
##moreThan1 21820
maxnum=0
maxindex=0
if(0):
    numOfMoreThan=0
    for line in open('..//data//bond.data'):
        bonds=line.split(',')
        ##print len(bonds)
        if len(bonds)>1:
            numOfMoreThan+=1
    print numOfMoreThan

numOFhuman=0
specis=dict()
specisNameP =r'OS=.*GN='


##r"(?<=<h1>).+?(?=<h1>)"
line='>sp|Q14524|SCN5A_HUMAN Sodium channel protein type 5 subunit alpha OS=Homo sapiens GN=SCN5A PE=1 SV=2'


for line in open('../data/sequene.fasta'):
    specisName = re.findall(specisNameP, line)
    if len(specisName)>0:
        specisName=specisName[0]
        specisName = re.findall(specisNameP, line)[0]
        specisName = specisName.replace("OS=", '')
        specisName = specisName.replace("GN=", '')
        specisName = specisName.replace(" ", '')
        if specisName in specis:
            specis[specisName]+=1
        else:
            specis[specisName]=1
specis=specis.items()
specis=[[v[1],v[0]] for v in specis]
specis.sort(reverse=True)

total=0
for value,key in specis:
    print(key,value)
    total=total+value
print()
print(total)


