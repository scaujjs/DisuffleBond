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
np.set_printoptions(threshold=np.nan)
import random
import os

from keras.models import Model
from keras.utils.vis_utils import plot_model
from keras import optimizers
np.set_printoptions(threshold=1000)



##6
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
if 0:
	modelName='modelForSStage_1'
	model=modelForSStage1()
	Xs = np.load(helper.pOOD + 'data/trainSStage1/X' + str(23) + '.data.npy')
	Ts=np.load(helper.pOOD + 'data/trainSStage1/T' + str(23) + '.data.npy')
	X_v =Xs
	T_v =Ts

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc=0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):

		model.fit_generator(helper.data_generaterForSStage1(128), epochs=100, steps_per_epoch=8*22 ,validation_data=[X_v,T_v],
							)##callbacks=[EarlyStopping(patience=3)]
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)

if 0:
	modelName = 'modelForSStage_1'
	model = modelForSStage1()
	Xs = np.load(helper.pOOD + 'data/trainSStage1/X' + str(24) + '.data.npy')
	Ts = np.load(helper.pOOD + 'data/trainSStage1/T' + str(24) + '.data.npy')
	X_v = Xs
	T_v = Ts

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	helper.testGet4metris(model,X_v,T_v)

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

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc=0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):

		model.fit_generator(helper.data_generaterForSStage1OneHot(128), epochs=100, steps_per_epoch=8*22 ,validation_data=[X_v,T_v],
							)##callbacks=[EarlyStopping(patience=3)]
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)

if 0:
	modelName = 'modelForSStageOnehot_1'
	model = modelForSStage1OneHot()
	Xs = np.load(helper.pOOD + 'data/trainSStage1OneHot/X' + str(24) + '.data.npy')
	Ts = np.load(helper.pOOD + 'data/trainSStage1OneHot/T' + str(24) + '.data.npy')
	X_v = Xs
	T_v = Ts

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	helper.testGet4metris(model,X_v,T_v)


def modelForSStage3():
	model = Sequential()
	model.add(Masking(mask_value=0, input_shape=(19, 100)))
	model.add(Bidirectional(LSTM(50,dropout=0.2,kernel_regularizer=regularizers.l2(0.005))))
	model.add((Dense(2, activation='softmax')))
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['acc'])
	return model


##this block used to triain
if 0:
	modelName='modelForSStage3'
	model=modelForSStage3()

	X_v = np.load(helper.pOOD + 'data/trainSStage3/Xv.data.npy')
	T_v=np.load(helper.pOOD + 'data/trainSStage3/Tv.data.npy')

	Xs=np.load(helper.pOOD + 'data/trainSStage3/X.data.npy')
	Ts=np.load(helper.pOOD + 'data/trainSStage3/T.data.npy')

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc=0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):
		model.fit(Xs,Ts,128, epochs=30,validation_data=[X_v,T_v]
							,callbacks=[EarlyStopping(patience=3)])
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)


if 0:
	modelName = 'modelForSStage3'
	model = modelForSStage3()
	Xs=np.load(helper.pOOD + 'data/trainSStage3/Xt.data.npy')
	Ts=np.load(helper.pOOD + 'data/trainSStage3/Tt.data.npy')
	X_v = Xs
	T_v = Ts

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	helper.testGet4metris(model,X_v,T_v)


if 0:
	modelName = 'modelForSStage1id30_word2vec'
	model = modelForSStage1()
	plot_model(model,"modelOfStage1",show_shapes=1)

	X = np.load(helper.pOOD + 'data/trainSStage1id30/X.data.npy')
	T = np.load(helper.pOOD + 'data/trainSStage1id30/T.data.npy')

	Xs=X[0:2500]
	Ts=T[0:2500]

	X_v=X[2500:3000]
	T_v=T[2500:3000]



	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):
		model.fit(Xs, Ts, 32, epochs=30, validation_data=[X_v, T_v]
				  , callbacks=[EarlyStopping(patience=3)])
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)

##model DSSP w2v
if 0:
	modelName = 'modelForSStage1DSSP_word2vec'
	model = modelForSStage1()
	##plot_model(model,"modelOfStage1",show_shapes=1)

	X = np.load(helper.pOOD + 'data/trainSStage1DSSP/X.data.npy')
	T = np.load(helper.pOOD + 'data/trainSStage1DSSP/T.data.npy')

	Xs=X[0:2500]
	Ts=T[0:2500]

	X_v=X[2500:3000]
	T_v=T[2500:3000]





	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	##model.load_weights(helper.pOOD + 'MODEL/keras/blstm_diff_ssb_v2d_3')
	for i in range(1):
		model.fit(Xs, Ts, 128, epochs=30, validation_data=[X_v, T_v]
				  , callbacks=[EarlyStopping(patience=3)])
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)


if 0:
	modelName = 'modelForSStage1DSSP_word2vec'
	model = modelForSStage1()
	##plot_model(model,"modelOfStage1",show_shapes=1)

	X = np.load(helper.pOOD + 'data/trainSStage1DSSP/X.data.npy')
	T = np.load(helper.pOOD + 'data/trainSStage1DSSP/T.data.npy')

	Xs = X[0:2500]
	Ts = T[0:2500]

	X_v = X[2500:3000]
	T_v = T[2500:3000]

	X_t = X[3000:3200]
	T_t = T[3000:3200]

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	helper.testGet4metris(model, X_v, T_v)

##model DSSP onehot
if 0:
	modelName = 'modelForSStage1DSSP_onehot'
	model = modelForSStage1OneHot()
	##plot_model(model,"modelOfStage1",show_shapes=1)

	X = np.load(helper.pOOD + 'data/trainSStage1OneHotDSSP/X.data.npy')
	T = np.load(helper.pOOD + 'data/trainSStage1OneHotDSSP/T.data.npy')

	Xs=X[0:2200]
	Ts=T[0:2200]

	X_v=X[2200:2800]
	T_v=T[2200:2800]

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName), callbacks=[EarlyStopping(patience=3)]
	maxacc = 0
	model.summary()
	for i in range(1):
		model.fit(Xs, Ts, 128, epochs=30, validation_data=[X_v, T_v])
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)

##model DSSP pssm
if 0:
	modelName = 'modelForSStage1DSSPpssm'
	model = modelForSStage1OneHot()
	##plot_model(model,"modelOfStage1",show_shapes=1)

	X = np.load(helper.pOOD + 'data/trainSStage1DSSPpssm/X.data.npy')
	T = np.load(helper.pOOD + 'data/trainSStage1DSSPpssm/T.data.npy')
	print(X)

	Xs=X[0:2200]
	Ts=T[0:2200]

	X_v=X[2200:2800]
	T_v=T[2200:2800]

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName), callbacks=[EarlyStopping(patience=3)]
	maxacc = 0
	model.summary()
	for i in range(1):
		model.fit(Xs, Ts, 128, epochs=30, validation_data=[X_v, T_v])
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)


def modelForpssmwin15():
## regular ,dropout=0.2,recurrent_regularizer=regularizers.l2(0.005),kernel_regularizer=regularizers.l2(0.005),
							 ##bias_regularizer=regularizers.l2(0.005),recurrent_dropout=0.2
	model = Sequential()
	model.add(Masking(mask_value=0, input_shape=(15, 20)))
	model.add(Bidirectional(LSTM(50,dropout=0.2,recurrent_regularizer=regularizers.l2(0.005),kernel_regularizer=regularizers.l2(0.005),
							 bias_regularizer=regularizers.l2(0.005),recurrent_dropout=0.2)))
	model.add((Dense(2, activation='softmax')))
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['acc'])
	return model


##model DSSP pssm win 15 &  culled,kernel_regularizer=regularizers.l2(0.01)
if 1:

    modelName = 'modelForpssmwindow151'
    model = modelForpssmwin15()
    Xs = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culled/X.data.npy')
    Ts = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culled/T.data.npy')

    X_v = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culled/Xt.data.npy')
    T_v = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culled/Tt.data.npy')



    ##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)##, callbacks=[EarlyStopping(patience=3)]
    maxacc = 0
    model.summary()
    for i in range(1):
        model.fit(Xs, Ts, 128, epochs=40, validation_data=[X_v, T_v],class_weight={0:1,}
                          )##, callbacks=[EarlyStopping(patience=3)])
        ##crightratio,rbondrightratio=test()
        model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True
                           )


if 0:
	modelName = 'modelForpssmwindow151'
	model = modelForpssmwin15()
	Xs=np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culled/Xt.data.npy')
	Ts=np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culled/Tt.data.npy')
	X_v=Xs
	T_v=Ts

	model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	helper.testGet4metris(model,X_v,T_v)
def modelForpssmwin15fn():

	model = Sequential()
	model.add(Dense(100, input_shape=(300,)))
	model.add(Dropout(0.2))
	model.add(Dense(2, activation='softmax'))
	model.compile(optimizer='rmsprop',
				  loss='categorical_crossentropy',
				  metrics=['accuracy'])
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['acc'])
	return model
##model DSSP pssm win 15 full connect
if 0:
    modelName = 'modelForpssmwindow15fullconnet'
    model = modelForpssmwin15fn()
    ##model=modelFullConnect15()
    ##plot_model(model,"modelOfStage1",show_shapes=1)

    X = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/X.data.npy')
    T = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/T.data.npy')
    Xs=X[0:25000]
    Ts=T[0:25000]
    Xs=Xs.reshape((Xs.shape[0],Xs.shape[1]*Xs.shape[2]))
    ##X_v = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/X_t.data.npy')
    ##T_v = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15culledBaseOnP/T_t.data.npy')
    X_v=X[25000:]
    T_v=T[25000:]
    X_v=X_v.reshape((X_v.shape[0],X_v.shape[1]*X_v.shape[2]))
	    ##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName), callbacks=[EarlyStopping(patience=3)]
    maxacc = 0
    model.summary()
    for i in range(1):
        model.fit(Xs, Ts, 128, epochs=100, validation_data=[X_v, T_v],callbacks=[EarlyStopping(patience=3)])
        ##crightratio,rbondrightratio=test() , callbacks=[EarlyStopping(patience=3)]
        model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)




def modelForpssmwin15fnRecurrent():

	model = Sequential()
	model.add(Masking(mask_value=0, input_shape=(50, 300)))
	model.add(Bidirectional(LSTM(50,return_sequences=True,dropout=0.2,recurrent_regularizer=regularizers.l2(0.005),kernel_regularizer=regularizers.l2(0.005),
								 bias_regularizer=regularizers.l2(0.005),recurrent_dropout=0.2)))
	model.add(TimeDistributed(Dense(2, activation='softmax')))
	model.compile(loss='categorical_crossentropy', optimizer='adam',metrics=['acc'])

	return model


##model DSSP pssm win 15 full connect+recurrent
if 0:
	modelName = 'modelForpssmwindow15fullconnetRecurrent'
	model = modelForpssmwin15fnRecurrent()
	##model=modelFullConnect15()
	##plot_model(model,"modelOfStage1",show_shapes=1)

	X = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15recurrent/X.data.npy')
	T = np.load(helper.pOOD + 'data/trainDSSPpssmwindow15recurrent/T.data.npy')
	##print(X)
	##X=X.reshape((X.shape[0],X.shape[1]*X.shape[2]))
	print(X.shape)
	Xs=X[0:2300]
	Ts=T[0:2300]

	X_v=X[2300:2700]
	T_v=T[2300:2700]

	##model.load_weights(helper.pOOD + 'MODEL/keras/' + modelName)
	maxacc = 0
	model.summary()
	for i in range(1):
		model.fit(Xs, Ts, 128, epochs=30, validation_data=[X_v, T_v]
				  , callbacks=[EarlyStopping(patience=3)])
		##crightratio,rbondrightratio=test()
		model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True)
