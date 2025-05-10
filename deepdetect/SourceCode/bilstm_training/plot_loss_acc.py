# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 11:29:39 2022

@author: hp
"""


from keras.callbacks import Callback
import matplotlib.pyplot as plt
import numpy as np


class LossHistoryOfBiLSTM(Callback):
    def on_train_begin(self, logs={}):
        self.losses = {'Batch': [], 'Epoch': []}
        self.accuracy = {'Batch': [], 'Epoch': []}
        self.val_loss = {'Batch': [], 'Epoch': []}
        self.val_acc = {'Batch': [], 'Epoch': []}
    
    def on_batch_end(self, batch, logs={}):
        self.losses['Batch'].append(logs.get('loss'))
        self.accuracy['Batch'].append(logs.get('binary_accuracy'))
        self.val_loss['Batch'].append(logs.get('val_loss'))
        self.val_acc['Batch'].append(logs.get('val_binary_accuracy'))
    
    def on_epoch_end(self, batch, logs={}):
        self.losses['Epoch'].append(logs.get('loss'))
        self.accuracy['Epoch'].append(logs.get('binary_accuracy'))
        self.val_loss['Epoch'].append(logs.get('val_loss'))
        self.val_acc['Epoch'].append(logs.get('val_binary_accuracy'))
    
    def loss_plot(self, loss_type, fig_path):
        iters = range(len(self.losses[loss_type]))
        # create a figure
        plt.figure(figsize=(8, 6), dpi=400)
        # lstm_acc
        plt.plot(iters, self.accuracy[loss_type], 'r--', 
                 label='train acc', linewidth=2)
        # loss
        plt.plot(iters, self.losses[loss_type], 'gray', 
                 label='train loss', linewidth=2)
        # validation
        if loss_type == 'Epoch':
            # val lstm acc
            plt.plot(iters, self.val_acc[loss_type], 'b--', 
                     label='val acc', linewidth=2)
            # val_loss
            plt.plot(iters, self.val_loss[loss_type], 'k', 
                     label='val loss', linewidth=2)
        # plot
        plt.grid(b=1, color='#CDC9C9', ls='-.')
        plt.xticks(fontsize=18)
        plt.ylim(0., 1.05)
        plt.yticks(np.linspace(0.,1.,6,endpoint=True), fontsize=18)
        plt.xlabel(loss_type, fontsize=20, labelpad=10)
        plt.ylabel('Loss/ACC', fontsize=20, labelpad=10)
        plt.legend(bbox_to_anchor=(0.72, 1.25), ncol=2, fontsize=18)
        plt.savefig(fig_path, dpi=400, bbox_inches='tight', pad_inches=0.2)
