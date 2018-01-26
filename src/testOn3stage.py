

##1 0 1 0 1 0 1 0 1 1 1 1 1 0 1 1 1 0 0 1 1 1 1 1 1 1 0 0 1 1 0 0 0 1 0 1 1
##1 1 0 0 1 0 1 0 0 1 1 1 1 0 1 1 0 0 1 0 0 1 1 0 1 1 0 1 1 0 0 0 0 0 0 0 1
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

trainExample=list()

index=0
tempfilename=''
tempseq=''
for line in open(helper.pOOD+'data/stage3Inp/Mix.txt'):
    if index%2==0:
        tempfilename=line.strip()
    else:
        tempseq=line.strip()
        trainExample.append((tempfilename,tempseq,[1,0]))
    index+=1

index=0
tempfilename=''
tempseq=''
for line in open(helper.pOOD+'data/stage3Inp/All.txt'):
    if index%2==0:
        tempfilename=line.strip()
    else:
        tempseq=line.strip()
        trainExample.append((tempfilename,tempseq,[0,1]))
    index+=1


random.shuffle(trainExample)
removelist=list()
for i in range(len(trainExample)):
    if len(trainExample[i][1])>600:
        removelist.append(trainExample[i])
for i in removelist:
    trainExample.remove(i)




def convert2word2vec(seq):


    input_length = len(seq) - 2
    X = np.ones((600, 100), np.float16) * (0)
    for i in range(input_length):
        aaw = helper.getaminoAcidword(seq, i + 2)

        X[i] = dictOfword2vect[aaw]

    return X

if 1:
    examples=list()
    for i in range(len(trainExample)):
        name=trainExample[i][0]
        seq=trainExample[i][1]
        Ti=trainExample[i][2]
        Xi=convert2word2vec(seq)
        examples.append((name,Xi,Ti))


    X_test=np.zeros((len(trainExample),600,100),np.float32)
    T_test=np.zeros((len(trainExample),2),np.float32)
    for i in range(len(trainExample)):
        seq=trainExample[i][1]
        Xi=convert2word2vec(seq)
        X_test[i]=Xi
        T_test[i]=trainExample[i][2]







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

    helper.testGet4metris(model,X_test,T_test)


def convert2Onthot(a1,a2,a3):
    X=np.zeros(66,np.float16)
    X[a1]=1
    X[a2+22]=1
    X[a3+44]=1
    return X
def convert2onehot(seq):



    input_length = len(seq) - 2
    X=np.zeros((600,66),np.float16)
    for i in range(input_length):
        aaw = helper.getaminoAcidword(seq, i + 2)
        a1=helper.aminoAcidDict.index(aaw[0])
        a2=helper.aminoAcidDict.index(aaw[1])
        a3=helper.aminoAcidDict.index(aaw[2])
        X[i]=convert2Onthot(a1, a2, a3)


    return X


if 1:
    X_test=np.zeros((261,600,66))
    T_test=np.zeros((261,2))


    for i in range(261):
        seq=trainExample[i][1]
        Xi=convert2onehot(seq)
        X_test[i]=Xi
        T_test[i]=trainExample[i][2]
        ##print(X_test[i],T_test[i])
    print(X_test.shape)
    print(T_test.shape)
    input_dim = 66
    output_dim = 2
    hidden_dim = 50
    batch_size = 128
    model_name = 'modelForSStageOnehot_1'

    model = Sequential()
    model.add(Masking(mask_value=0, input_shape=(600, 66)))
    model.add(Bidirectional(LSTM(50)))
    model.add((Dense(output_dim, activation='softmax')))
    model.compile(loss='categorical_crossentropy', optimizer='adam',
                  metrics=['acc'])
    model.load_weights(helper.pOOD + 'MODEL/keras/' + model_name)

    helper.testGet4metris(model,X_test,T_test)

