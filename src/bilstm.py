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

##model DSSP onehot code, and do the stage 1 prediction---- whether its a full proteins or mix proteins
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

##model DSSP pssm, this block also try to pick the mix from full, but the code method is PSSM
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



## this model use BiLSTM, and context of 15 aa the predict the state of CYS in the centre of the
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

## this block is used train, the model used here is @modelForpssmwin15
## this is a stange thing which confused my team, the acc of validation based on CYS suffle is higher(92)
## that the acc based on proteins suffle
if 0:

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
        model.fit(Xs, Ts, 128, epochs=40, validation_data=[X_v, T_v],class_weight={0:1,1:3}
                          )##, callbacks=[EarlyStopping(patience=3)])
        ##crightratio,rbondrightratio=test()
        model.save_weights(helper.pOOD + 'MODEL/keras/' + modelName, overwrite=True
                           )


## this block is used to test,with the model @modelForpssmwin15
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


## this model used tranditional full-connect
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
