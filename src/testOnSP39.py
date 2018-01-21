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
            print(error)
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
    input_length = len(seq) - 2
    X = np.ones((600, 100), np.float16) * (0)
    T = np.ones((600,3), np.float16) * (0)
    T[input_length-1:]=[0,0,0]
    for i in range(input_length):
        aaw = helper.getaminoAcidword(seq, i + 2)

        X[i] = dictOfword2vect[aaw]

        if seq[i+1]=='C':
            T[i+1]=[0,0,1]
        else:
            T[i+1]=[1,0,0]

    for position1,position2 in pairs:
        T[position1-1]=[0,1,0]
        T[position2-1]=[0,1,0]

    ##print("------------")
    ##for i in range(len(seq)):
        ##print(T[i].argmax(),seq[i])
        ##np.savetxt(pOOD+'data/train/temp_'+str(num)+'x.data',X,fmt='%10.5f',delimiter=',')
        ##np.savetxt(pOOD+'data/train/temp_'+str(num)+'t.data',T,fmt='%10.5f',delimiter=',')
    return X, T



X_test=np.zeros((446,600,100))
T_test=np.zeros((446,600,3))


for i in range(417):
    seq=dictOfProtein[listOfkey[i]]
    bond=dictOfBond[listOfkey[i]]
    Xi,Ti=convert2word2vec(seq,bond)
    X_test[i]=Xi
    T_test[i]=Ti


input_dim = 100
output_dim = 3
hidden_dim = 50
batch_size = 128
model_name = 'blstm_diff_ssb_v2d_3'

model = Sequential()
model.add(Masking(mask_value=0, input_shape=(600, 100)))
model.add(Bidirectional(LSTM(50, return_sequences=True)))
model.add(TimeDistributed(Dense(output_dim, activation='softmax')))
model.compile(loss='categorical_crossentropy', optimizer='adam',
              metrics=['acc', helper.single_class_accuracyC, helper.single_class_accuracy])
model.load_weights(helper.pOOD + 'MODEL/keras/' + model_name)

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