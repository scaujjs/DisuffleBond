
import os

import PIL.Image as image
import numpy as np
import random
np.set_printoptions(threshold=np.NaN)
import tensorflow as tf
import math
pOOD=os.getcwd()+'/../'
##215768(this is the num of Cys totally)
##143796(this is the num of Cys in bond)
## 71972(this is the num of Cys not in bond)

##protein num=25582
##protein all Cys in Bond=8258
##      50988Cys in bond; 0 Cys notInBond
##protein mix Cys 17324
##      92808Cys in bond; 71972 Cys notInBond
##the different type of amino acid
##['M', 'K', 'V', 'L', 'I', 'S', 'P', 'N', 'G', 'A', 'R', 'H', 'D', 'Q', 'W', 'E', 'C', 'F', 'T', 'Y', 'X', 'U']
##there 22 type of amino acid



##in the Dataset from DSSP-unique
##full is 2017
##mix is 1182


##total Cys is 16736
##13040 is bonded
##3692 is free




##third Dataset from PIESES
## 13018 raw
##filt CYS<2 7666
##pdb file 7425
##dssp file 7420
## get 7654 seq finally
## after filt odd bond 7479 left

## after filt length more 1000
## 7427 left

##total CYS:35325 27655
##bonded CYS:7018 6996



aminoAcidDict=['M', 'K', 'V', 'L', 'I', 'S', 'P', 'N', 'G', 'A', 'R', 'H', 'D', 'Q', 'W', 'E', 'C', 'F', 'T', 'Y', 'X', 'U']

##this block use to red protein word2vec
dictOfword2vect=dict()
num=0
for line in open(pOOD+'data/protVec_100d_3grams.csv'):
    vector=np.zeros(100,float)
    items=line.split('\t')
    for j in range(100):
        vector[j]=float(items[j+1].replace('\"',''))
    dictOfword2vect[items[0].replace('\"','')]=vector


##this
def parseLineOfBond(line,reduceReplicate=False,reduceLine=False):
    flagOfPass=False
    items=line.split(',')
    key=items[0]
    pairs=list()
    items.remove(items[0])
    for item in items:
        stringPair=item.split('_')
        pair=(int(stringPair[0]),int(stringPair[1]))
        if reduceReplicate and pair[0]==pair[1]:
            flagOfPass=True
            pass
        else:
            pairs.append(pair)
    if flagOfPass and reduceLine:
        return (key,[])
    else:
        return (key,pairs)

##in this block I fillter seq longer than 600, stanger bond, interchain bond
## last 24000 seq
def generateNewRawData():
    dictOfSequence=dict()
    dictOfPair=dict()
    num=0
    def detectError(seq):
        errorFlag=0
        for i in range(len(seq)-2):
            aaw = getaminoAcidword(seq, i+2)
            if aaw not in dictOfword2vect:
                #print aaw
                errorFlag =1
                break
        return errorFlag
    for line in open(pOOD+'NewData/newSeqFromPDB_unique'):

        items=line.split()[0].split(',')
        if len(items[1])>599:
            continue
        if detectError(items[1]):
            continue
        if len(items[1])<50:
            continue
        num+=1
        dictOfSequence[items[0]]=items[1]

    num=0
    numOfPairs=0
    for line in open(pOOD+'NewData/newBondsFromPDB_unique'):
        (key,pairs)=parseLineOfBond(line,reduceReplicate=True,reduceLine=True)
        if len(pairs)>0:
            dictOfPair[key]=pairs
            numOfPairs+=len(pairs)
            num+=1


    for key,_ in dictOfSequence.items():
        if key not in dictOfPair:
            dictOfSequence.pop(key)
    for key,_ in dictOfPair.items():
        if key not in dictOfSequence:
            dictOfPair.pop(key)



    listOfSeqAndBond=list()
    listOfP=dictOfPair.items()
    for i in range(len(listOfP)):
        key=listOfP[i][0]

        lineOfSeq=key+','+dictOfSequence[key]+'\n'
        lineOfPair=key+','


        pairs=dictOfPair[key]
        for j in range(len(pairs)):
            lineOfPair+=str(pairs[j][0])+'_'+str(pairs[j][1])+','
        lineOfPair=lineOfPair[0:-1]+'\n'

        listOfSeqAndBond.append((lineOfSeq,lineOfPair))
    random.shuffle(listOfSeqAndBond)
    random.shuffle(listOfSeqAndBond)
    random.shuffle(listOfSeqAndBond)

    listOfnewSeq=list()
    listOfnewPair=list()


    for lineOfnewSeq,lineOfnewPair in listOfSeqAndBond:
        listOfnewSeq.append(lineOfnewSeq)
        listOfnewPair.append(lineOfnewPair)


    fileOfNewSeq=open(pOOD+'data/seqDSSP.data','w')
    fileOfNewSeq.writelines(listOfnewSeq)
    fileOfNewPair=open(pOOD+'data/bondDSSP.data','w')
    fileOfNewPair.writelines(listOfnewPair)

##
def findAllCword(seq,length):
    aAWs=dict()
    for i in range(len(seq)):
        if seq[i]=='C':
            aaw=getaminoAcidword(seq,i+1,length)
            if aaw in aAWs:
                aAWs[aaw]+=1
            else:
                aAWs[aaw]=1

    return aAWs


## the index here is the index in fasta, which means 1 is the frist rather than 0
def getaminoAcidword(seq,index,length=3):
    aminoAcidword=''
    start=index-1-length//2
    for i in range(length):
        if 0<=start+i<=len(seq)-1:
            aminoAcidword+=seq[start+i]
        else:
            aminoAcidword+='X'
    return aminoAcidword


def loadTrainData():
    listOfKey=list()
    listOfSeq=list()
    listOfPairs=list()

    for line in open(pOOD+'data/newSeq.data'):
        key,seq=line.strip().split(',')
        listOfSeq.append(seq)
        listOfKey.append(key)

    for line in open(pOOD+'data/newBond.data'):
        items=line.strip().split(',')
        pairs=list()
        for i in range(len(items)-1):
            strpair=items[i+1].split('_')
            pair=(int(strpair[0]),int(strpair[1]))
            pairs.append(pair)
        listOfPairs.append(pairs)

    return(listOfKey,listOfSeq,listOfPairs)


def loadTrainDataId30():
    listOfKey=list()
    listOfSeq=list()
    listOfPairs=list()

    for line in open(pOOD+'data/seqId30.data'):
        key,seq=line.strip().split(',')
        listOfSeq.append(seq)
        listOfKey.append(key)

    for line in open(pOOD+'data/bondId30.data'):
        items=line.strip().split(',')
        pairs=list()
        for i in range(len(items)-1):
            strpair=items[i+1].split('_')
            pair=(int(strpair[0]),int(strpair[1]))
            pairs.append(pair)
        listOfPairs.append(pairs)

    return(listOfKey,listOfSeq,listOfPairs)

def loadTrainDataDSSP():
    listOfKey=list()
    listOfSeq=list()
    listOfPairs=list()

    for line in open(pOOD+'data/seqDSSP.data'):
        key,seq=line.strip().split(',')
        listOfSeq.append(seq)
        listOfKey.append(key)

    for line in open(pOOD+'data/bondDSSP.data'):
        items=line.strip().split(',')
        pairs=list()
        for i in range(len(items)-1):
            strpair=items[i+1].split('_')
            pair=(int(strpair[0]),int(strpair[1]))
            pairs.append(pair)
        listOfPairs.append(pairs)

    return(listOfKey,listOfSeq,listOfPairs)

def loadTrainDataDSSPculled25():
    listOfKey=list()
    listOfSeq=list()
    listOfPairs=list()

    for line in open(pOOD+'data/newSeqFromPDBculled25'):
        key,seq=line.strip().split(',')
        listOfSeq.append(seq)
        listOfKey.append(key)

    for line in open(pOOD+'data/newBondsFromPDBculled25'):
        items=line.strip().split(',')
        pairs=list()
        for i in range(len(items)-1):
            strpair=items[i+1].split('_')
            pair=(int(strpair[0]),int(strpair[1]))
            pairs.append(pair)
        listOfPairs.append(pairs)

    return(listOfKey,listOfSeq,listOfPairs)

def loadTrainDataDSSPculled25Dict():
    listOfSeq=dict()
    listOfPairs=dict()

    for line in open(pOOD+'data/newSeqFromPDBculled25'):
        key,seq=line.strip().split(',')
        listOfSeq[key]=seq


    for line in open(pOOD+'data/newBondsFromPDBculled25'):
        items=line.strip().split(',')
        pairs=list()
        for i in range(len(items)-1):
            strpair=items[i+1].split('_')
            pair=(int(strpair[0]),int(strpair[1]))
            pairs.append(pair)

        listOfPairs[items[0]]=pairs

    return(listOfSeq,listOfPairs)


##this task predict the odd there is a ssBond
##[ nornml,inpair,notinpair,empy]
def generateTrainingDataSSbondOdd():
    def convert2word2vec(seq, pairs):
        input_length = len(seq) - 2
        X = np.ones((600, 100), np.float16) * (0)
        T = np.ones((600,3), np.float16) * (0)
        T[input_length-1:]=[0,0,0]
        for i in range(input_length):
            aaw = getaminoAcidword(seq, i + 2)

            X[i] = dictOfword2vect[aaw]

            if seq[i+1]=='C':
                T[i+1]=[0,0,1]
            else:
                T[i+1]=[1,0,0]

        for position1,position2 in pairs:
            T[position1-1]=[0,1,0]
            T[position2-1]=[0,1,0]

            ##np.savetxt(pOOD+'data/train/temp_'+str(num)+'x.data',X,fmt='%10.5f',delimiter=',')
            ##np.savetxt(pOOD+'data/train/temp_'+str(num)+'t.data',T,fmt='%10.5f',delimiter=',')
        return X, T

    (listOfKey, listOfSeq, listOfPairs)=loadTrainData()

    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs_Ts=np.zeros((1024,600,103),np.float16)




    index=0


    for i in range(len(listOfPairs)):


        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        (X,T)=convert2word2vec(seq,pairs)


        Xs_Ts[(index)%1024,:,0:100]=X[:,:]
        Xs_Ts[(index)%1024,:,100:]=T[:,:]
        index+=1

        if (index) % 1024 == 0:
            print('num ' + str(i))
            np.save(pOOD+'data/trainSSOdd/'+str(index//1024)+'.data',Xs_Ts)

def data_generatorForSsOdd(batch_size):
    files=os.listdir(pOOD+'data/trainSSOdd/')
    indexlist1024=list()
    indexlist24=list()
    for i in range(1024):
        indexlist1024.append(i)
    for i in range(23):
        indexlist24.append(i)
    X=np.zeros((batch_size,600,100),np.float)
    T=np.zeros((batch_size,600,3),np.float)
    random.shuffle(indexlist24)

    fileNum=26
    while True:
        for i in range(22):
            ##Xs_Ts=np.load(pOOD+'data/trainSSOdd/'+str(24)+'.data.npy')
            Xs_Ts=np.load(pOOD+'data/trainSSOdd/'+str(indexlist24[i]+1)+'.data.npy')

            random.shuffle(indexlist1024)
            for j in range(1024//batch_size):
                for k in range(batch_size):
                    X[k]=Xs_Ts[indexlist1024[j*batch_size+k],:,0:100]
                    T[k]=Xs_Ts[indexlist1024[j*batch_size+k],:,100:]

                yield X,T

def generteTrainingDataSStage1():
    def convert2word2vec(seq, pairs):
        tempCys=0
        for am in seq:
            if am=='C':
                tempCys+=1
        if tempCys==len(pairs)*2:
            T=1
        else:
            T=0


        input_length = len(seq) - 2
        X = np.ones((600, 100), np.float16) * (0)
        for i in range(input_length):
            aaw = getaminoAcidword(seq, i + 2)

            X[i] = dictOfword2vect[aaw]

        return X, T

    (listOfKey, listOfSeq, listOfPairs)=loadTrainData()
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs=np.zeros((1024,600,100),np.float16)
    Ts=np.zeros((1024,2),np.float16)
    index=0

    for i in range(len(listOfPairs)):
        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        (X,T)=convert2word2vec(seq,pairs)

        Xs[(index)%1024,:,:]=X[:,:]
        if T==0:

            Ts[(index)%1024]=[1,0]
        else:
            Ts[(index)%1024]=[0,1]
        index+=1

        if (index) % 1024 == 0:
            print('num ' + str(i))
            np.save(pOOD+'data/trainSStage1/X'+str(index//1024)+'.data',Xs)
            np.save(pOOD+'data/trainSStage1/T'+str(index//1024)+'.data',Ts)

def data_generaterForSStage1(batch_size):
    indexlist1024=list()
    indexlist22=list()
    for i in range(1024):
        indexlist1024.append(i)
    for i in range(22):
        indexlist22.append(i)
    X=np.zeros((batch_size,600,100),np.float)
    T=np.zeros((batch_size,2),np.float)
    random.shuffle(indexlist22)
    print(T.shape)

    while True:
        for i in range(22):
            ##Xs_Ts=np.load(pOOD+'data/trainSSOdd/'+str(24)+'.data.npy')
            Xs=np.load(pOOD+'data/trainSStage1/X'+str(indexlist22[i]+1)+'.data.npy')
            Ts=np.load(pOOD+'data/trainSStage1/T'+str(indexlist22[i]+1)+'.data.npy')
            random.shuffle(indexlist1024)
            for j in range(1024//batch_size):
                for k in range(batch_size):
                    X[k]=Xs[indexlist1024[j*batch_size+k],:,:]
                    T[k]=Ts[indexlist1024[j*batch_size+k],:]


                yield X,T
                                          
def generateTrainingDataSSStage1OneHot():
    def convert2Onthot(a1,a2,a3):
        X=np.zeros(66,np.float16)
        X[a1]=1
        X[a2+22]=1
        X[a3+44]=1
        return X
    def convert2word2vec(seq, pairs):
        tempCys=0
        for am in seq:
            if am=='C':
                tempCys+=1
        if tempCys==len(pairs)*2:
            T=1
        else:
            T=0


        input_length = len(seq) - 2
        X=np.zeros((600,66),np.float16)
        for i in range(input_length):
            aaw = getaminoAcidword(seq, i + 2)
            a1=aminoAcidDict.index(aaw[0])
            a2=aminoAcidDict.index(aaw[1])
            a3=aminoAcidDict.index(aaw[2])
            X[i]=convert2Onthot(a1, a2, a3)


        return X, T

    (listOfKey, listOfSeq, listOfPairs)=loadTrainData()
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs=np.zeros((1024,600,66),np.float16)
    Ts=np.zeros((1024,2),np.float16)
    index=0

    for i in range(len(listOfPairs)):
        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        (X,T)=convert2word2vec(seq,pairs)

        Xs[(index)%1024,:,:]=X[:,:]
        if T==0:

            Ts[(index)%1024]=[1,0]
        else:
            Ts[(index)%1024]=[0,1]
        index+=1

        if (index) % 1024 == 0:
            print('num ' + str(i))
            np.save(pOOD+'data/trainSStage1OneHot/X'+str(index//1024)+'.data',Xs)
            np.save(pOOD+'data/trainSStage1OneHot/T'+str(index//1024)+'.data',Ts)

def data_generaterForSStage1OneHot(batch_size):

    indexlist1024=list()
    indexlist22=list()
    for i in range(1024):
        indexlist1024.append(i)
    for i in range(22):
        indexlist22.append(i)
    X=np.zeros((batch_size,600,66),np.float)
    T=np.zeros((batch_size,2),np.float)
    random.shuffle(indexlist22)
    print(T.shape)

    while True:
        for i in range(22):
            ##Xs_Ts=np.load(pOOD+'data/trainSSOdd/'+str(24)+'.data.npy')
            Xs=np.load(pOOD+'data/trainSStage1OneHot/X'+str(indexlist22[i]+1)+'.data.npy')
            Ts=np.load(pOOD+'data/trainSStage1OneHot/T'+str(indexlist22[i]+1)+'.data.npy')
            random.shuffle(indexlist1024)
            for j in range(1024//batch_size):
                for k in range(batch_size):
                    X[k]=Xs[indexlist1024[j*batch_size+k],:,:]
                    T[k]=Ts[indexlist1024[j*batch_size+k],:]


                yield X,T

def generateTrainingDataSSStage3():
    ##the index here is fasta format, which mean 1 is the first letter
    ## and the value of the window size must be odd
    def getSubSeq(seq,index,windowSize=21):
        subSeq=''
        haldWindowsize=windowSize//2
        for i in range(windowSize):
            tempIndex=index-1-haldWindowsize+i
            if 0<=tempIndex<=len(seq)-1:
                subSeq=subSeq+seq[tempIndex]
            else:
                subSeq+='X'
        return subSeq

    def getExampleInASeq(seq, pairs):
        example=list()
        tempCys=0
        for am in seq:
            if am=='C':
                tempCys+=1
        if tempCys==len(pairs)*2:
            pass
        else:
            for i in range(len(seq)):
                if seq[i]=='C':
                    index=i+1
                    T=[1,0]
                    for a,b in pairs:
                        if index==a or index==b:
                            T=[0,1]
                            break
                    subSeq=getSubSeq(seq,index)
                    example.append((subSeq,T))
        return example

    def convertWord2Vec(example):
        finaalyexample=list()
        #print('_____________')
        for seq,t in example:
            ##print(seq,t)
            errorFlag=0
            X = np.ones((19, 100), np.float16) * (0)
            for i in range(len(seq)-2):
                aaw = getaminoAcidword(seq, i + 2)
                if aaw not in dictOfword2vect:
                    errorFlag=1
                    break
                X[i] = dictOfword2vect[aaw]
            if not errorFlag:
                finaalyexample.append((X,t))
        random.shuffle(finaalyexample)
        return finaalyexample

    (listOfKey, listOfSeq, listOfPairs)=loadTrainData()
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs=np.zeros((100000,19,100),np.float16)
    Ts=np.zeros((100000,2),np.float16)
    Xs_v=np.zeros((30000,19,100),np.float16)
    Ts_v=np.zeros((30000,2),np.float16)
    Xs_t=np.zeros((27000,19,100),np.float16)
    Ts_t=np.zeros((27000,2),np.float16)

    index=0

    for i in range(len(listOfPairs)):
        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        example=getExampleInASeq(seq,pairs)
        finalExm=convertWord2Vec(example)
        for X,T in finalExm:
            ##print(X,T)
            if index<100000:
                Xs[index,:,:]=X
                Ts[index]=T
            if 100000<=index<130000:
                Xs_v[index-100000,:,:]=X
                Ts_v[index-100000]=T
            if 130000<=index<157000:
                Xs_t[index-130000,:,:]=X
                Ts_t[index-130000]=T

            index += 1

            if index%10000==0:
                print(index)
    np.save(pOOD+'data/trainSStage3/X.data',Xs)
    np.save(pOOD+'data/trainSStage3/T.data',Ts)

    np.save(pOOD+'data/trainSStage3/Xv.data',Xs_v)
    np.save(pOOD+'data/trainSStage3/Tv.data',Ts_v)

    np.save(pOOD+'data/trainSStage3/Xt.data',Xs_t)
    np.save(pOOD+'data/trainSStage3/Tt.data',Ts_t)



def generteTrainingDataSStage1id30():
    def convert2word2vec(seq, pairs):
        tempCys=0
        for am in seq:
            if am=='C':
                tempCys+=1
        if tempCys==len(pairs)*2:
            T=1
        else:
            T=0


        input_length = len(seq) - 2
        X = np.ones((600, 100), np.float16) * (0)
        for i in range(input_length):
            aaw = getaminoAcidword(seq, i + 2)

            X[i] = dictOfword2vect[aaw]

        return X, T

    (listOfKey, listOfSeq, listOfPairs)=loadTrainDataId30()
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs=np.zeros((3474,600,100),np.float16)
    Ts=np.zeros((3474,2),np.float16)
    index=0

    for i in range(len(listOfPairs)):
        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        (X,T)=convert2word2vec(seq,pairs)

        Xs[index,:,:]=X[:,:]
        if T==0:

            Ts[index]=[1,0]
        else:
            Ts[index]=[0,1]
        index+=1


    np.save(pOOD+'data/trainSStage1id30/X'+'.data',Xs)
    np.save(pOOD+'data/trainSStage1id30/T'+'.data',Ts)

def generteTrainingDataSStage1DSSP():
    def convert2word2vec(seq, pairs):
        tempCys=0
        for am in seq:
            if am=='C':
                tempCys+=1
        if tempCys==len(pairs)*2:
            T=1
        else:
            T=0


        input_length = len(seq) - 2
        X = np.ones((600, 100), np.float16) * (0)
        for i in range(input_length):
            aaw = getaminoAcidword(seq, i + 2)

            X[i] = dictOfword2vect[aaw]

        return X, T

    (listOfKey, listOfSeq, listOfPairs)=loadTrainDataDSSP()
    print(len(listOfKey))
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs=np.zeros((len(listOfKey),600,100),np.float16)
    Ts=np.zeros((len(listOfKey),2),np.float16)
    index=0

    for i in range(len(listOfPairs)):
        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        (X,T)=convert2word2vec(seq,pairs)

        Xs[index,:,:]=X[:,:]
        if T==0:

            Ts[index]=[1,0]
        else:
            Ts[index]=[0,1]
        index+=1

    np.save(pOOD+'data/trainSStage1DSSP/X'+'.data',Xs)
    np.save(pOOD+'data/trainSStage1DSSP/T'+'.data',Ts)

def generateTraininDataStage1DSSPpssm():
    def convert2word2vec(seq, pairs):
        tempCys=0
        for am in seq:
            if am=='C':
                tempCys+=1
        if tempCys==len(pairs)*2:
            T=1
        else:
            T=0
        return T
    Xraw=np.load(pOOD + 'data/seqDSSPpssm.npy')
    listOfKey,listOfSeq,listOfPairs=loadTrainDataDSSP()
    Xs=np.zeros((len(listOfKey),600,60),np.float16)
    Ts=np.zeros((len(listOfKey),2),np.float16)

    def convert20To60(X):
        Xnew=np.zeros((600,60),np.float32)
        one=X[0]
        two=X[1]
        for i in range(598):
            temp=X[i+2]
            if temp.max()==0 and temp.min()==0:
                break
            else:
                Xnew[i,0:20]=one
                Xnew[i,20:40]=two
                Xnew[i,40:]=temp
                one=two
                two=temp
        return Xnew



    for i in range(len(listOfPairs)):
        seq=listOfSeq[i]
        pairs=listOfPairs[i]
        T=convert2word2vec(seq,pairs)
        X=Xraw[i]
        X=convert20To60(X)
        if T==0:
            Ts[i]=[1,0]
        else:
            Ts[i]=[0,1]
        Xs[i]=X
    print(Ts)

    print(Xs.shape,Ts.shape)
    np.save(pOOD+'data/trainSStage1DSSPpssm/X'+'.data',Xs)
    np.save(pOOD+'data/trainSStage1DSSPpssm/T'+'.data',Ts)

##generateTraininDataStage1DSSPpssm()

def generateTrainingDataSSStage1OneHotDSSP():
    def convert2Onthot(a1, a2, a3):
        X = np.zeros(66, np.float16)
        X[a1] = 1
        X[a2 + 22] = 1
        X[a3 + 44] = 1
        return X

    def convert2word2vec(seq, pairs):
        tempCys = 0
        for am in seq:
            if am == 'C':
                tempCys += 1
        if tempCys == len(pairs) * 2:
            T = 1
        else:
            T = 0

        input_length = len(seq) - 2
        X = np.zeros((600, 66), np.float16)
        for i in range(input_length):
            aaw = getaminoAcidword(seq, i + 2)
            a1 = aminoAcidDict.index(aaw[0])
            a2 = aminoAcidDict.index(aaw[1])
            a3 = aminoAcidDict.index(aaw[2])
            X[i] = convert2Onthot(a1, a2, a3)

        return X, T

    (listOfKey, listOfSeq, listOfPairs) = loadTrainDataDSSP()
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not
    Xs = np.zeros((len(listOfKey), 600, 66), np.float16)
    Ts = np.zeros((len(listOfKey), 2), np.float16)
    index = 0

    for i in range(len(listOfPairs)):
        seq = listOfSeq[i]
        pairs = listOfPairs[i]
        (X, T) = convert2word2vec(seq, pairs)

        Xs[index, :, :] = X[:, :]
        if T == 0:

            Ts[index] = [1, 0]
        else:
            Ts[index] = [0, 1]
        index += 1


    np.save(pOOD + 'data/trainSStage1OneHotDSSP/X.data', Xs)
    np.save(pOOD + 'data/trainSStage1OneHotDSSP/T.data', Ts)


def generatrTrainingData15pssm():
    ##conf ## lengthOfSubseq
    lengthOfSubseq = 15
    halfwin = lengthOfSubseq // 2

    def generateSubexample(seq, pssm, index, pairs):
        X = np.zeros((lengthOfSubseq, 20), np.float32)

        for i in range(lengthOfSubseq):
            tempindex = index - halfwin + i
            if tempindex < 0 or tempindex > len(seq) - 1:
                X[i] = np.zeros(20, np.float32)
            else:
                X[i] = pssm[i]

        T = [1., 0.]
        for p1, p2 in pairs:
            if index == p1 - 1 or index == p2 - 1:
                T = [0., 1.]
                break
        return X, T

    Xraw = np.load(pOOD + 'data/seqDSSPpssm.npy')
    listOfKey, listOfSeq, listOfPairs = loadTrainDataDSSP()

    examples = list()
    for i in range(len(listOfKey)):
        seq = listOfSeq[i]
        for j in range(len(seq)):
            if seq[j] == 'C':
                X, T = generateSubexample(seq, Xraw[i], j, listOfPairs[i])
                print(T)
                examples.append((X, T))
    random.shuffle(examples)
    Xs = np.zeros((len(examples), lengthOfSubseq, 20), np.float32)
    Ts = np.zeros((len(examples), 2), np.float32)

    numOf01 = 0
    numOf10 = 0
    for i in range(len(examples)):
        if examples[i][1] == [0., 1.]:
            numOf01 += 1
        if examples[i][1] == [1., 0.]:
            numOf10 += 1
        Xs[i] = examples[i][0]
        Ts[i] = examples[i][1]

    print(numOf01)
    print(numOf10)
    print(len(examples))

    np.save(pOOD + 'data/trainDSSPpssmwindow15/X.data', Xs)
    np.save(pOOD + 'data/trainDSSPpssmwindow15/T.data', Ts)

def generatrTrainingData15pssmculled():
    ##conf ## lengthOfSubseq
    lengthOfSubseq = 15
    halfwin = lengthOfSubseq // 2

    def generateSubexample(seq, pssm, index, pairs):
        X = np.zeros((lengthOfSubseq, 20), np.float32)

        for i in range(lengthOfSubseq):
            tempindex = index - halfwin + i
            if tempindex < 0 or tempindex > len(seq) - 1:
                X[i] = np.zeros(20, np.float32)
            else:
                X[i] = pssm[i]

        T = [1., 0.]
        for p1, p2 in pairs:
            if index == p1 - 1 or index == p2 - 1:
                T = [0., 1.]
                break
        return X, T

    Xraw = np.load(pOOD + 'data/seqDSSPpssmculled25.npy')
    listOfSeq, listOfPairs = loadTrainDataDSSPculled25Dict()
    ##print(listOfPairs.items())
    listOfIdentity=list()
    for line in open(pOOD + 'data/indexOfCulled25'):
        listOfIdentity.append(line.strip())

    examples = list()
    for i in range(len(listOfIdentity)):
        seq = listOfSeq[listOfIdentity[i]]
        for j in range(len(seq)):
            if seq[j] == 'C':
                X, T = generateSubexample(seq, Xraw[i], j, listOfPairs[listOfIdentity[i]])
                print(T)
                ##print(X)
                examples.append((X, T))

    random.shuffle(examples)
    Xs = np.zeros((len(examples), lengthOfSubseq, 20), np.float32)
    Ts = np.zeros((len(examples), 2), np.float32)

    numOf01 = 0
    numOf10 = 0
    for i in range(len(examples)):
        if examples[i][1] == [0., 1.]:
            numOf01 += 1
        if examples[i][1] == [1., 0.]:
            numOf10 += 1
        Xs[i] = examples[i][0]

        Ts[i] = examples[i][1]

    print(numOf01)
    print(numOf10)
    print(len(examples))

    np.save(pOOD + 'data/trainDSSPpssmwindow15culled/X.data', Xs)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culled/T.data', Ts)
##generatrTrainingData15pssmculled()

def generatrTrainingData15pssmculledBasedOnProteins():
    ##conf ## lengthOfSubseq
    lengthOfSubseq = 15
    halfwin = lengthOfSubseq // 2

    def generateSubexample(seq, pssm, index, pairs):
        X = np.zeros((lengthOfSubseq, 20), np.float32)

        for i in range(lengthOfSubseq):
            tempindex = index - halfwin + i
            if tempindex < 0 or tempindex > len(seq) - 1:
                X[i] = np.zeros(20, np.float32)
            else:
                X[i] = pssm[i]

        T = [1., 0.]
        for p1, p2 in pairs:
            if index == p1 - 1 or index == p2 - 1:
                T = [0., 1.]
                break
        return X, T

    Xraw = np.load(pOOD + 'data/seqDSSPpssmculled25.npy')
    listOfSeq, listOfPairs = loadTrainDataDSSPculled25Dict()
    ##print(listOfPairs.items())
    listOfIdentity=list()
    for line in open(pOOD + 'data/indexOfCulled25'):
        listOfIdentity.append(line.strip())

    ##random.shuffle(listOfIdentity)
    onefifth=len(listOfIdentity)//5


    examplesOfTrain = list()
    for i in range(onefifth*4):
        seq = listOfSeq[listOfIdentity[i]]
        pair=listOfPairs[listOfIdentity[i]]
        for j in range(len(seq)):
            if seq[j] == 'C':
                X, T = generateSubexample(seq, Xraw[i], j, pair)
                examplesOfTrain.append((X, T))
    random.shuffle(examplesOfTrain)
    Xs = np.zeros((len(examplesOfTrain), lengthOfSubseq, 20), np.float32)
    Ts = np.zeros((len(examplesOfTrain), 2), np.float32)
    numOf01 = 0
    numOf10 = 0
    for i in range(len(examplesOfTrain)):
        if examplesOfTrain[i][1] == [0., 1.]:
            numOf01 += 1
        if examplesOfTrain[i][1] == [1., 0.]:
            numOf10 += 1
        Xs[i] = examplesOfTrain[i][0]
        Ts[i] = examplesOfTrain[i][1]

    print(numOf01)
    print(numOf10)
    print(len(examplesOfTrain))



    listOfTestkey=list()
    examplesOfTrst = list()
    for i in range(len(listOfIdentity)-onefifth*4):
        seq = listOfSeq[listOfIdentity[i+onefifth*4]]
        pair=listOfPairs[listOfIdentity[i+onefifth*4]]
        listOfTestkey.append(listOfIdentity[i+onefifth*4]+'\n')
        for j in range(len(seq)):
            if seq[j] == 'C':
                X, T = generateSubexample(seq, Xraw[i+onefifth*4], j, pair)
                ##print(T)
                ##print(X)
                examplesOfTrst.append((X, T))
    random.shuffle(examplesOfTrst)
    X_t = np.zeros((len(examplesOfTrst), lengthOfSubseq, 20), np.float32)
    T_t = np.zeros((len(examplesOfTrst), 2), np.float32)

    numOf01 = 0
    numOf10 = 0
    for i in range(len(examplesOfTrst)):
        if examplesOfTrst[i][1] == [0., 1.]:
            numOf01 += 1
        if examplesOfTrst[i][1] == [1., 0.]:
            numOf10 += 1
        X_t[i] = examplesOfTrst[i][0]
        print(examplesOfTrst[i][0])
        T_t[i] = examplesOfTrst[i][1]

    print(numOf01)
    print(numOf10)
    print(len(examplesOfTrst))

    f=open(pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/idOftest','w')
    f.writelines(listOfTestkey)
    f.close()
    np.save(pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/X_t.data', X_t)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/T_t.data', T_t)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/X.data', Xs)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/T.data', Ts)
##generatrTrainingData15pssmculledBasedOnProteins()

def generatrTrainingData15pssmculledBasedOnProteinsMethod2():
    ##conf ## lengthOfSubseq
    lengthOfSubseq = 15
    halfwin = lengthOfSubseq // 2

    def generateSubexample(seq, pssm, index, pairs):
        X = np.zeros((lengthOfSubseq, 20), np.float32)

        for i in range(lengthOfSubseq):
            tempindex = index - halfwin + i
            if tempindex < 0 or tempindex > len(seq) - 1:
                X[i] = np.zeros(20, np.float32)
            else:
                X[i] = pssm[i]

        T = [1., 0.]
        for p1, p2 in pairs:
            if index == p1 - 1 or index == p2 - 1:
                T = [0., 1.]
                break
        return X, T

    Xraw = np.load(pOOD + 'data/seqDSSPpssmculled25.npy')
    listOfSeq, listOfPairs = loadTrainDataDSSPculled25Dict()
    ##print(listOfPairs.items())
    listOfIdentity = list()
    for line in open(pOOD + 'data/indexOfCulled25'):
        listOfIdentity.append(line.strip())

    examples = list()
    for i in range(len(listOfIdentity)):
        seq = listOfSeq[listOfIdentity[i]]
        for j in range(len(seq)):
            if seq[j] == 'C':
                X, T = generateSubexample(seq, Xraw[i], j, listOfPairs[listOfIdentity[i]])
                print(T)
                ##print(X)
                examples.append((X, T))
    TRAINNUM=len(examples)//5*4
    TestNUM=len(examples)-TRAINNUM
    testExample=examples[0:TestNUM]

    trainExample=examples[TestNUM:]

    random.shuffle(trainExample)
    random.shuffle(testExample)

    Xs = np.zeros((len(trainExample), lengthOfSubseq, 20), np.float32)
    Ts = np.zeros((len(trainExample), 2), np.float32)
    numOf01 = 0
    numOf10 = 0
    for i in range(len(trainExample)):
        if trainExample[i][1] == [0., 1.]:
            numOf01 += 1
        if trainExample[i][1] == [1., 0.]:
            numOf10 += 1
        Xs[i] = trainExample[i][0]

        Ts[i] = trainExample[i][1]
    print(numOf01)
    print(numOf10)
    print(len(trainExample))
    Xt = np.zeros((len(testExample), lengthOfSubseq, 20), np.float32)
    Tt = np.zeros((len(testExample), 2), np.float32)
    numOf01 = 0
    numOf10 = 0
    for i in range(len(testExample)):
        if testExample[i][1] == [0., 1.]:
            numOf01 += 1
        if testExample[i][1] == [1., 0.]:
            numOf10 += 1
        Xt[i] = testExample[i][0]

        Tt[i] = testExample[i][1]
    print(numOf01)
    print(numOf10)
    print(len(testExample))

    np.save(pOOD + 'data/trainDSSPpssmwindow15culled/X.data', Xs)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culled/T.data', Ts)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culled/Xt.data', Xt)
    np.save(pOOD + 'data/trainDSSPpssmwindow15culled/Tt.data', Tt)
##generatrTrainingData15pssmculledBasedOnProteinsMethod2()
def generatrTrainingData15pssmRecurrent():

##conf ## lengthOfSubseq
    lengthOfSubseq=15
    halfwin=lengthOfSubseq//2
    def generateSubexample(seq,pssm,index,pairs):
        X=np.zeros((lengthOfSubseq,20),np.float32)

        for i in range(lengthOfSubseq):
            tempindex=index-halfwin+i
            if tempindex<0 or tempindex>len(seq)-1:
                X[i]=np.zeros(20,np.float32)
            else:
                X[i]=pssm[i]

        T=[1.,0.]
        for p1,p2 in pairs:
            if index==p1-1 or index==p2 -1:
                T=[0.,1.]
                break
        X=X.reshape(X.shape[0]*X.shape[1])

        return X,T




    Xraw=np.load(pOOD + 'data/seqDSSPpssm.npy')
    listOfKey,listOfSeq,listOfPairs=loadTrainDataDSSP()

    Xs=np.zeros((2824,50,300),np.float32)
    Ts=np.zeros((2825,50,2),np.float32)
    for i in range(len(listOfKey)):
        XC=np.zeros((50,300),np.float32)
        TC=np.zeros((50,2),np.float32)
        seq=listOfSeq[i]
        index=0
        for j in range(len(seq)):
            if seq[j]=='C':
                X,T=generateSubexample(seq,Xraw[i],j,listOfPairs[i])
                XC[index]=X
                TC[index]=T
                index+=1
        Xs[i]=XC
        Ts[i]=TC
        print(TC)



    np.save(pOOD + 'data/trainDSSPpssmwindow15recurrent/X.data', Xs)
    np.save(pOOD + 'data/trainDSSPpssmwindow15recurrent/T.data', Ts)

##generatrTrainingData15pssm()
##generatrTrainingData15pssmRecurrent()

from keras import backend as K

INTERESTING_CLASS_ID = 1  # Choose the class of interest

def single_class_accuracy(y_true, y_pred):
    INTERESTING_CLASS_ID1 = 1  # Choose the class of interest
    INTERESTING_CLASS_ID2=2
    y_true=y_true*tf.constant([0,1,1],dtype=tf.float32)
    y_pred=y_pred*tf.constant([0,1,1],dtype=tf.float32)
    class_id_true = K.argmax(y_true, axis=2)


    class_id_true=tf.reshape(class_id_true,[76800])
    class_id_preds = K.argmax(y_pred, axis=2)
    class_id_preds=tf.reshape(class_id_preds,[76800])
    # Replace class_id_preds with class_id_true for recall here
    accuracy_mask1 = K.cast(K.equal(class_id_true, INTERESTING_CLASS_ID1), dtype=tf.float32)
    class_acc_tensor1 = K.cast(K.equal(class_id_true, class_id_preds),dtype=tf.float32) * accuracy_mask1

    accuracy_mask2 = K.cast(K.equal(class_id_true, INTERESTING_CLASS_ID2), dtype=tf.float32)
    class_acc_tensor2 = K.cast(K.equal(class_id_true, class_id_preds),dtype=tf.float32) * accuracy_mask2

    acc = (K.sum(class_acc_tensor1)+K.sum(class_acc_tensor2)) / (K.maximum(K.sum(accuracy_mask1), 1)+K.maximum(K.sum(accuracy_mask2), 1))



    return acc


def single_class_accuracyC(y_true, y_pred):
    INTERESTING_CLASS_ID1 = 1  # Choose the class of interest
    class_id_true = K.argmax(y_true, axis=2)


    class_id_true=tf.reshape(class_id_true,[76800])
    class_id_preds = K.argmax(y_pred, axis=2)
    class_id_preds=tf.reshape(class_id_preds,[76800])
    # Replace class_id_preds with class_id_true for recall here
    accuracy_mask1 = K.cast(K.equal(class_id_true, INTERESTING_CLASS_ID1), dtype=tf.float32)
    class_acc_tensor1 = K.cast(K.equal(class_id_true, class_id_preds),dtype=tf.float32) * accuracy_mask1

    acc = (K.sum(class_acc_tensor1)) / (K.maximum(K.sum(accuracy_mask1), 1))



    return ac


if 0:
    a=np.array([[[0, 1, 0], [0, 0.4, 0.3],[0, 0.4, 0.3]],
                [[ 0.1, 2, 1], [0,  1, 2],[0, 0.4, 0.3]]])

    y_true=tf.constant(a)

    y_pred=tf.constant([[[0, 1, 0], [0, 0, 0.3],[0, 0.4, 0.3]],
                        [[ 0.1, 1, 1], [ 3, 1, 2],[0, 0.4, 0.3]],
                         ])



    with tf.Session() as sess:
        class_id_true = K.argmax(y_true, axis=2)
        class_id_true = tf.reshape(class_id_true, [6])
        print(class_id_true.eval())
        class_id_preds = K.argmax(y_pred, axis=2)
        class_id_preds= tf.reshape(class_id_preds, [6])
        print(class_id_preds.eval())

        result=single_class_accuracy(y_true, y_pred)
        print(result.eval())
        accuracy_mask = K.cast(K.equal(class_id_true, INTERESTING_CLASS_ID), dtype=tf.float32)
        print(accuracy_mask.eval())

##statistic description
if 0 :
    if 0:
        numOfProtein=0
        for line in open(pOOD+'data/newBond.data'):
            numOfProtein+=1
        print(numOfProtein)
        numOfProtein=0
        for line in open(pOOD + 'data/newSeq.data'):
            numOfProtein+=1
        print(numOfProtein)

    if 0:
        Cysinbond = 0

        for line in open(pOOD+'data/newBond.data'):
            item=line.split(',')
            Cysinbond+=(len(item)-1)*2

        print(Cysinbond)
    if 0:

        CysTotal=0
        for line in open(pOOD + 'data/newSeq.data'):
            item = line.split(',')
            for am in item[1]:
                if am=="C":
                    CysTotal+=1

        print(CysTotal)
    if 0:
        CysInFull=0
        listOfKey,listOfSeq,listOfPair=loadTrainData()
        fullCysInBond=0
        for i in range(len(listOfSeq)):

            tempOFTotalCys=0
            tempofBondingCys=2*len(listOfPair[i])
            seq=listOfSeq[i]
            for am in seq:
                if am=='C':
                    tempOFTotalCys+=1
            if tempOFTotalCys==tempofBondingCys:
                fullCysInBond+=1
                CysInFull+=tempOFTotalCys
        print(fullCysInBond)
        print(CysInFull)

    if 1:
        typeAm=list()
        for line in open(pOOD + 'data/newSeq.data'):
            item = line.split(',')
            for am in item[1]:
                if am not in typeAm:
                    typeAm.append(am)
        print(typeAm)
        print(len(typeAm))


def testGet4metris(model,X_v,T_v):
    T_v=np.argmax(T_v,1)
    Y_v = model.predict(X_v)
    print(Y_v.shape)
    Y_v = np.argmax(Y_v, 1)
    print(Y_v.shape)
    print(T_v.shape)
    print(T_v)
    print(Y_v)

    TP = 0
    FP = 0
    FN = 0
    TN = 0

    for i in range(X_v.shape[0]):
        if T_v[i] == 0 and Y_v[i] == 0:
            TP += 1
        if T_v[i] == 1 and Y_v[i] == 1:
            TN += 1
        if T_v[i] == 1 and Y_v[i] == 0:
            FP += 1
        if T_v[i] == 0 and Y_v[i] == 1:
            FN += 1
    print(TP,FP,FN,TN)
    acc = (TP + TN) * 1.0 / X_v.shape[0]
    sensitivity = TP * 1.0 / (TP + FN)
    specificity = TN * 1.0 / (TN + FP)
    MCC = (TP * TN - FP * FN) * 1.0 / pow((TP + FN) * 1.0 * (TP + FP) * (TN + FP) * (TN + FN), 0.5)
    print(acc, sensitivity, specificity, MCC)
    return (acc, sensitivity, specificity, MCC)

##generateNewRawData()
##generteTrainingDataSStage1id30()

##generateTrainingDataSSStage1OneHotDSSP()
##generteTrainingDataSStage1DSSP()
if 0:
    T = np.load(pOOD + 'data/trainSStage1OneHotDSSP/T.data.npy')
    T=np.argmax(T,1)
    numOfmix=0
    numOfFull=0
    for i in range(T.shape[0]):
        if T[i]==0:
            numOfmix+=1
        if T[i]==1:
            numOfFull+=1
    print(numOfmix,numOfFull)


if 0:
    listOfKey, listOfSeq, listOfPairs = loadTrainDataDSSP()

    totalCys=0
    for seq in listOfSeq:
        for aa in seq:
            if aa=='C':
                totalCys+=1
    bondCys=0
    for pairs in listOfPairs:
        bondCys+=2*(len(pairs))

    print(totalCys,bondCys)


if 0:
    import os
    import time

    while 1:
        time.sleep(5)
        for i in range(82):
            stri=str(i)
            name=''
            if len(stri)==1:
                name='0'+stri
            else:
                name=stri
            filename='nr.'+name+'.tar.gz'
            if os.path.exists(filename):
                if os.path.exists(filename+'.aria2'):
                    print("start recover"+filename)
                    command='aria2c -x 4 -j 6 -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/'+filename
                    os.system(command)
                else:
                    pass
            else:
                print("start "+filename)
                command = 'aria2c -x 4 -j 6 -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/' + filename
                os.system(command)


if 0:
    listOfKey,listOfSeq,listOfPairs=loadTrainDataDSSPculled25()

    for i in range(len(listOfKey)):
        ##print(i)
        listOfshowBefore=list()
        pathhh=pOOD+'fasta2pssm/fasta2/'+listOfKey[i]+'.fasta'
        if pathhh not  in listOfshowBefore:
            listOfshowBefore.append(pathhh)
        else:
            print(pathhh)
        if 0:
            fw=open(pOOD+'fasta2pssm/fasta2/'+listOfKey[i]+'.fasta','w')
            fw.write('>'+listOfKey[i]+'\n')
            fw.write(listOfSeq[i]+'\n')
            fw.close()

## THIS bLOCK is used to use psiblast
if 0:
    files=os.listdir(pOOD+'fasta2pssm/fasta2/')
    index=0
    for file in files:
        pssmname=file[0:-5]+'pssm'
        inputp=pOOD+'fasta2pssm/fasta2/'+file
        outputp=pOOD+'fasta2pssm/pssm2/'+pssmname
        command='psiblast -inclusion_ethresh 0.001 -num_iterations 3 -num_threads 5 -db swissprot -in_msa '+inputp+' -out_ascii_pssm '+outputp
        os.system(command)
        index+=1
        print(file+' have been generated: '+str(index))
 

def convert2standard(pssmRaw):
    maxV = pssmRaw.max()
    minV = pssmRaw.min()
    meanV = pssmRaw.mean()
    errorflag=0
    if maxV == minV:
        return np.zeros(20,np.float32),0
    else:
        for i in range(20):
            pssmRaw[i] = (pssmRaw[i] - meanV) / (maxV - minV)
            if math.isnan(pssmRaw[i]):
                errorflag=1
        return pssmRaw,errorflag


def parePSSM2matrix(filePath):
    start = 0
    pattern = 'A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V'
    PSSMlist = list()
    for line in open(filePath):

        if start:
            pssm20column = np.zeros(20, np.float32)
            usefulitem = list()
            items = line.split(" ")
            for item in items:
                if item.strip() != '':
                    usefulitem.append(item.strip())
            if len(usefulitem) != 44:
                break
            for i in range(20):
                num=int(usefulitem[i + 2])
                if num>1000 or num<-1000:
                    print(num)
                    print(filePath)
                    print('error')
                pssm20column[i]=num
            PSSMlist.append(pssm20column)
        if pattern in line:
            start = 1
            continue
    if len(PSSMlist)>1000:
        return None
    else:
        PSSM = np.zeros((1000, 20), np.float32)
        for i in range(len(PSSMlist)):
            ##standize PSSM
            standardPSSM,errorFlag = convert2standard(PSSMlist[i])
            if errorFlag:
                print('error')
                print(i)
                print(filePath)
                break

            PSSM[i] = standardPSSM
    ##print(len(PSSMlist))
    return PSSM


##parse PSSM
def generatePSSMraw():

    listOfSeq,listOfPairs=loadTrainDataDSSPculled25Dict()
    files=list()


    for key,seq in listOfSeq.items():
        files.append(key+'.pssm')
    print(files)
    listPSSMraw=list()
    listOfIdentity=list()
    index=0
    for file in files:
        ##print(file)
        filePath=pOOD+'fasta2pssm/pssm2/'+file
        PSSM=parePSSM2matrix(filePath)
        if PSSM is not None:
            listPSSMraw.append(PSSM)

            listOfIdentity.append(file[:-5]+'\n')
            index+=1
            print(file[-11:-5])
            print(index)

        else:
            print('return none')



    X=np.zeros((len(listOfIdentity),1000,20),np.float32)
    for i in range(len(listOfIdentity)):
        X[i]=listPSSMraw[i]

    f=open(pOOD + 'data/indexOfCulled25','w')
    f.writelines(listOfIdentity)
    np.save(pOOD + 'data/seqDSSPpssmculled25', X)




##generatePSSMraw()
if 0:
    listOfKey,listOfSeq,listOfPairs=loadTrainDataDSSP()
    maxC=0
    A=1
    biggerThanA=0
    for seq in listOfSeq:
        numinthisS=0
        for aa in seq:
            if aa=="C":
                numinthisS+=1
        if numinthisS>maxC:
            maxC=numinthisS
        if numinthisS> A:
            biggerThanA+=1
    print(maxC)
    print(biggerThanA)


##this block generate query 1 to search agaist PDB
if 0:
    print('hello world')
    seqDict=dict()
    start=0
    pidTemp=''
    seqTemp=''
    lineNum=100
    for line in open(pOOD+'data/PIESES/fasta.txt'):
        if line.startswith('>'):
            if start:
                seqDict[pidTemp]=seqTemp
                ##print(pidTemp)
                ##print(seqTemp)
                seqTemp=''
            else:
                start=1
            pidTemp=line.strip().replace("|PDBID|CHAIN|SEQUENCE",'').replace(">",'').lower()
        else:
            seqTemp+=line.strip()
    seqDict[pidTemp]=seqTemp
    import os

    idList = list()
    for line in open(pOOD+'data/PIESES/pid'):
        if line.startswith('IDs'):
            continue

        items = line.split(' ')
        for item in items:
            if item.strip() != '':
                idList.append(item.lower())
                ##print(item.lower())
                break
    for pid,_ in seqDict.items():

        if pid.replace(':','') not in idList:
            ##print(pid.replace(':',''))
            seqDict.pop(pid)

    print(len(seqDict.items()))
    print(len(idList))

    for pid,seq in seqDict.items():
        numOfCYS=0
        for aa in seq:
            if aa=='C':
                numOfCYS+=1
        if numOfCYS<2:
            seqDict.pop(pid)
    print(len(seqDict.items()))

    idShowBefore=list()
    idWithChainId=list()
    for key,seq in seqDict.items():
        pid=key[0:4]
        if pid not in idShowBefore:
            idShowBefore.append(pid)
    print(len(idShowBefore))

    for key,seq in seqDict.items():
        idWithChainId.append(key+'\n')

    f= open('../listOfCulled','w')
    f.writelines(idWithChainId)
    f.close()



    if 0:
        query = list()
        index = 0
        strtemp = ''
        for pid, seq in seqDict.items():
            strtemp += (pid[0:4] + ',')
            index += 1
            if index % 20 == 0 and index > 1:
                strtemp += '\n'
                query.append(strtemp)
                strtemp = ''
        query.append(strtemp)
        f = open('query1', 'w')
        f.writelines(query)
        f.close()



## this block count seq less or equal than 1000 is  7432
if 0:
    listOfKey=list()
    listOfSeq=list()
    listOfPairs=list()
    numOflessThan1000=0
    for line in open(pOOD+'data/newSeqFromPDBculled25'):
        key,seq=line.strip().split(',')
        listOfSeq.append(seq)
        listOfKey.append(key)
        if len(seq)<=1000:
            numOflessThan1000+=1
    print(numOflessThan1000)

## this block count how many CYS, how many Bonded CYS
if 0:
    totalCys=0
    bondedCys=0
    listOfKey,listOfSeq,listOfPairs=loadTrainDataDSSPculled25()
    for seq in listOfSeq:
        for aa in seq:
            if aa=='C':
                totalCys+=1
    for pairs in listOfPairs:
        bondedCys+=2*(len(pairs))

    print(totalCys)
    print(bondedCys)


