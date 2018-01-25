import helper
import numpy as np
from keras.models import Model
from keras.models import Sequential
from keras.layers import TimeDistributed
from keras.layers import Bidirectional
from keras.layers import LSTM,Masking
from keras.callbacks import EarlyStopping
from keras import regularizers
from keras.layers import Input,Conv2D, MaxPooling2D,UpSampling2D,Dense, Dropout
from keras.layers.core import  Activation,  Flatten, Reshape
import random

dictOfword2vect=helper.dictOfword2vect
dictOfBond=dict()
dictOfProtein=dict()
listOfkey=list()


i=-1
error=0
name=''
for line in open(helper.pOOD+'test/SP39.fasta'):
    i+=1
    if i%2==0:
        name=line.strip().replace('>','')
    else:
        seq=line.strip()
        if len(seq)>600:
            error+=1
            ##print(error)
            continue
        listOfkey.append(name)
        dictOfProtein[name]=seq




for key,_ in dictOfProtein.items():
    bond=list()
    i=0
    for line in open(helper.pOOD+'test/SP39_BOND/'+key):
        if i==0:
            i+=1
            continue
        item=line.split(' ')
        bond.append((int(item[0]),int(item[1])))
    dictOfBond[key]=bond

def convert2word2vec(seq, pairs):
    numOfC=0
    for am in seq:
        if am=="C":
            numOfC+=1
    if numOfC==len(pairs)*2:
        T=[0,1]
    else:
        T=[1,0]


    input_length = len(seq) - 2
    X = np.ones((600, 100), np.float16) * (0)
    for i in range(input_length):
        aaw = helper.getaminoAcidword(seq, i + 2)

        X[i] = dictOfword2vect[aaw]


    ##print("------------")
    ##for i in range(len(seq)):
        ##print(T[i].argmax(),seq[i])
        ##np.savetxt(pOOD+'data/train/temp_'+str(num)+'x.data',X,fmt='%10.5f',delimiter=',')
        ##np.savetxt(pOOD+'data/train/temp_'+str(num)+'t.data',T,fmt='%10.5f',delimiter=',')
    return X, T


if 1:
    X_test=np.zeros((417,600,100))
    T_test=np.zeros((417,2))


    for i in range(417):
        seq=dictOfProtein[listOfkey[i]]
        bond=dictOfBond[listOfkey[i]]
        Xi,Ti=convert2word2vec(seq,bond)
        X_test[i]=Xi
        T_test[i]=Ti

    print(X_test.shape)
    print(T_test.shape)
    input_dim = 100
    output_dim = 2
    hidden_dim = 50
    batch_size = 128
    model_name = 'modelForSStage_1'

    model = Sequential()
    model.add(Masking(mask_value=0, input_shape=(600, 100)))
    model.add(Bidirectional(LSTM(50)))
    model.add((Dense(output_dim, activation='softmax')))
    model.compile(loss='categorical_crossentropy', optimizer='adam',
                  metrics=['acc'])
    model.load_weights(helper.pOOD + 'MODEL/keras/' + model_name)
    if 0:
        totalCB = 0
        totalCNB = 0
        rihgtCB = 0
        rihgtCNB = 0
        P = model.predict(X_test)*[0,1,1]
        P = P.argmax(axis=-1)
        T = T_test.argmax(axis=-1)

        for i in range(417):
            print("_________")
            for j in range(T.shape[1]):
                if T[i][j]==0:
                    continue
                print(P[i][j],T[i][j])
        allrightNum = 0

        for i in range(417):
            flag = 1

            for j in range(600):
                if T[i][j] == 1:
                    totalCB += 1
                    if P[i][j] == 1:
                        rihgtCB += 1
                    else:
                        flag=0

                if T[i][j] == 2:
                    totalCNB += 1
                    if P[i][j] == 2:
                        rihgtCNB += 1
                    else:
                        flag=0
            if flag:
                allrightNum+=1

        print(totalCB, totalCNB, rihgtCB, rihgtCNB)

        crightratio = 1.0*(rihgtCB + rihgtCNB) / (totalCB + totalCNB)
        rbondrightratio = (1.0*rihgtCB / totalCB)
        print(crightratio)
        print(rbondrightratio)
        print(allrightNum/417.0)

    helper.testGet4metris(model,X_test,T_test)

if 0:
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
                aaw =helper.getaminoAcidword(seq, i + 2)
                if aaw not in dictOfword2vect:
                    errorFlag=1
                    break
                X[i] = dictOfword2vect[aaw]
            if not errorFlag:
                finaalyexample.append((X,t))
        random.shuffle(finaalyexample)
        return finaalyexample

    (listOfKey, listOfSeq, listOfPairs)=helper.loadTrainData()
    ## T is in the 1 of 101, 1 means there is a ss, 0 means not

    index=0
    numOfExam=0
    X=np.zeros((1693,19,100),np.float32)
    T=np.zeros((1693,2),np.float32)
    for i in range(417):
        seq=dictOfProtein[listOfkey[i]]
        bond=dictOfBond[listOfkey[i]]
        example=getExampleInASeq(seq,bond)
        finalExm=convertWord2Vec(example)
        numOfExam+=len(finalExm)
        for exm in  finalExm:
            X[index]=exm[0]
            T[index]=exm[1]
            index+=1
    model_name='modelForSStage3'
    model = Sequential()
    model.add(Masking(mask_value=0, input_shape=(19, 100)))
    model.add(Bidirectional(LSTM(50, kernel_regularizer=regularizers.l2(0.001))))
    model.add((Dense(2, activation='softmax')))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['acc'])
    model.load_weights(helper.pOOD + 'MODEL/keras/' + model_name)
    helper.testGet4metris(model,X,T)