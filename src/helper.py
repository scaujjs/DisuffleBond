
import os

import numpy as np
import random
np.set_printoptions(threshold=np.NaN)
import tensorflow as tf
import math
pOOD=os.getcwd()+'/../'
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


##  this method will convert the PSSM directory to 2 files. the first file is a npy file which store
## all standard pssm sequence and the second file is the index(PID) of pssm sequence
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


