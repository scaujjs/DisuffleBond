from keras.models import Model
from keras.models import Sequential
from keras.layers import TimeDistributed
from keras.layers import Bidirectional
from keras.layers import LSTM,Masking
from keras.callbacks import EarlyStopping
from keras import regularizers
from keras.layers import Input,Conv2D, MaxPooling2D,UpSampling2D,Dense, Dropout
from keras.layers.core import  Activation,  Flatten, Reshape
import helper
import numpy as np
import random
import os
from keras.models import Model

from keras import optimizers
np.set_printoptions(threshold=1000)



## shangjie has been here
input_dim = 100
output_dim = 3
hidden_dim = 50
batch_size = 128
model_name = "blstm_diff_ssb_v2d"

##create the model
def modelForSSBOnd():

	model = Sequential()
	model.add(Masking(mask_value=0, input_shape=(600, 100)))
	model.add(Bidirectional(LSTM(50, return_sequences=True)))
	model.add(TimeDistributed(Dense(output_dim, activation='softmax')))
	model.compile(loss='categorical_crossentropy', optimizer='adam',metrics=[helper.single_class_accuracy])
	return model


##this is the test function, this function will be called when test the model
def test():
	totalCB=0
	totalCNB=0
	rihgtCB=0
	rihgtCNB=0
	P=model.predict(X_test)
	P=P*[0,1,1]
	P=P.argmax(axis=-1)
	print(P)
	T=T_test.argmax(axis=-1)

	for i in range(1024):
		print("_________")
		for j in range(T.shape[1]):
			if T[i][j] == 0:
				continue
			print(P[i][j], T[i][j])

	allrightNum=0
	for i in range(1024):
		flagOfAllright=1

		for j in range(600):
			if T[i][j]==1:
				totalCB+=1
				if P[i][j]==1:
					rihgtCB+=1
				else:
					flagOfAllright=0

			if T[i][j]==2:
				totalCNB+=1
				if P[i][j]==2:
					rihgtCNB+=1
				else:
					flagOfAllright=0
		if flagOfAllright:
			allrightNum+=1


	crightratio=1.0*(rihgtCB+rihgtCNB)/(totalCB+totalCNB)
	rbondrightratio=1.0*rihgtCB/totalCB
	print(totalCB,totalCNB,rihgtCB,rihgtCNB)
	print(crightratio)
	##print(rbondrightratio)
	print((allrightNum+0.0)/1024)

	return (crightratio,rbondrightratio)


##this block used to test
if 0:
	model=modelForSSBOnd()
	model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')

	Xs_Ts = np.load(helper.pOOD + 'data/trainSSOdd/22.data.npy')
	X_test=Xs_Ts[:,:,0:100]
	T_test=Xs_Ts[:,:,100:]
	sequence=list()
	keys=list()
	index=1
	test()



##this block used to triain
if 0:
	model=modelForSSBOnd()
	Xs_Ts = np.load(helper.pOOD + 'data/trainSSOdd/' + str(23) + '.data.npy')
	X_v =Xs_Ts[:,:,0:100]
	T_v =Xs_Ts[:,:,100:]

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + model_name)
	maxacc=0
	model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):

		model.fit_generator(helper.data_generatorForSsOdd(128), epochs=70, steps_per_epoch=8 * 22,validation_data=[X_v,T_v],
							)##callbacks=[EarlyStopping(patience=3)]
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + model_name+"_3", overwrite=True)


def modelForSStage1():
	model = Sequential()
	model.add(Masking(mask_value=0, input_shape=(600, 100)))
	model.add(Bidirectional(LSTM(50)))
	model.add((Dense(2, activation='softmax')))
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['acc'])
	return model


##this block used to triain
if 1:
	modelName='modelForSStage_1'
	model=modelForSStage1()
	Xs = np.load(helper.pOOD + 'data/trainSStage1/X' + str(23) + '.data.npy')
	Ts=np.load(helper.pOOD + 'data/trainSStage1/T' + str(23) + '.data.npy')
	X_v =Xs
	T_v =Ts

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc=0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):

		model.fit_generator(helper.data_generaterForSStage1(128), epochs=10, steps_per_epoch=8*22 ,validation_data=[X_v,T_v],
							)##callbacks=[EarlyStopping(patience=3)]
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)


def modelForSStage1OneHot():
	model = Sequential()
	model.add(Masking(mask_value=0, input_shape=(600, 66)))
	model.add(Bidirectional(LSTM(50)))
	model.add((Dense(2, activation='softmax')))
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['acc'])
	return model


##this block used to triain
if 0:
	modelName='modelForSStageOnehot_1'
	model=modelForSStage1OneHot()
	Xs = np.load(helper.pOOD + 'data/trainSStage1OneHot/X' + str(23) + '.data.npy')
	Ts=np.load(helper.pOOD + 'data/trainSStage1OneHot/T' + str(23) + '.data.npy')
	X_v =Xs
	T_v =Ts

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc=0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):

		model.fit_generator(helper.data_generaterForSStage1OneHot(128), epochs=10, steps_per_epoch=8*22 ,validation_data=[X_v,T_v],
							)##callbacks=[EarlyStopping(patience=3)]
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)
