# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:51:18 2019

@author: hp
"""


import random as rn
import numpy as np
import tensorflow as tf
# import tensorflow.compat.v1 as tf
import os
from data_coding import data_encoding
from plot_loss_acc import LossHistoryOfBiLSTM
from keras.layers import (Input, Embedding, LSTM, Bidirectional,
                          normalization, Dense)
from keras.models import Model
#from keras.optimizers import adam
#from keras.utils import plot_model
import time


os.environ['PYTHONHASHSEED'] = str(0)
rn.seed(0)  # set python random seed (0)
np.random.seed(1)  # set numpy random seed (1)
# tf.disable_v2_behavior()
tf.get_logger().setLevel('ERROR')  # ignore warning, print error only
# set tensorflow random seed (2)
tf.set_random_seed(2)  # tf == 1.10.0
# tf.compat.v1.set_random_seed(2)  # tf == 1.15.0
# tf.random.set_seed(2)  # tf >= 2.0
# force to use single thread in tensorflow, considering that the reproduction
# can be effected by multithreading
# tf == 1.10.0
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, 
                              inter_op_parallelism_threads=1)
# # tf == 1.15.0
# session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1, 
#                                         inter_op_parallelism_threads=1)
# # tf >= 2.0
# session = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1, 
#                                     inter_op_parallelism_threads=1)
# session_conf = tf.compat.v1.Session(config=config)
os.environ['TF_DETERMINISTIC_OPS'] = '1'
os.environ['TF_CUDNN_DETERMINISTIC'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
# from keras.backend.tensorflow_backend import set_session
# config = tf.ConfigProto()
# config.gpu_options.per_process_gpu_memory_fraction = 0.3
# #config.gpu_options.allow_growth = True
# set_session(tf.Session(config=config))


# build model
def bilstm(max_len):
    print("Building model.")
    # model architecture
    seq_input = Input(shape=(max_len,),
                      name='Sequence')
    x = Embedding(input_dim=21,
                  output_dim=10,
                  input_length=max_len,
                  mask_zero=True)(seq_input)
    lstm_features = Bidirectional(LSTM(units=20,
                                       activation='softsign',
                                       dropout=0.5,
                                       unroll=True))(x)
    x = normalization.BatchNormalization()(lstm_features)
    lstm_predict = Dense(units=1,
                         activation='sigmoid',
                         name='final_predict')(x)
    model = Model(inputs=seq_input, outputs=lstm_predict)
    
    # model compile
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=['binary_accuracy'])
    
    # plot_model(model,
    #            to_file='BiLSTM_%s.png' % protease,
    #            show_shapes=True,
    #            show_layer_names=True)
    # model.summary()
    
    return model


# training
def training(pep_path, protease):
    # load training data
    x_train, y_train, x_val, y_val, max_len = data_encoding(pep_path=pep_path)
    print("There are %s peptides in the training data set."
          % (len(y_train) + len(y_val)))
    
    # build model
    model = bilstm(max_len)
    
    # create a LossHistory
    history = LossHistoryOfBiLSTM()
    
    # training and validation
    print("Training.")
    start = time.clock()
    model.fit(x_train, y_train,
              batch_size=64, epochs=50,
              verbose=2, callbacks=[history],
              validation_data=(x_val, y_val),
              shuffle=False, workers=1)
    end = time.clock()
    print("Time cost of training is %ss." % (end - start))
    
    # plot loss-acc figure
    history.loss_plot('Epoch', r'BiLSTM_LossACC_%s.png' % protease)
    
    # save model and weights
    print("Saving model and weights.")
    model_json = model.to_json()
    with open(r'BiLSTM_%s.json' % protease, 'w') as md:
        md.write(model_json)
    model.save_weights(r'BiLSTM_%s.h5' % protease)
    print("Model and weights have been saved.")


if __name__ == '__main__':
    trdata = '2019Xu_LysargiNase_Yeast'
    pep_path = '../datasets/%s.txt' % trdata
    protease = trdata.split('_')[1]
    print("| Training BiLSTM on %s |" % trdata)
    training(pep_path, protease)