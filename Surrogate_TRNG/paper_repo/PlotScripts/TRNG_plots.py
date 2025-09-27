#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 22:13:43 2025

@author: kevin

"""

# Generate the kramers time data and further generate the histograms

# Required modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import curve_fit
import scipy as sc
import configparser as conparse
import json
import pandas as pd

#%% Required Functions

# Theoretical plot functions
def F(x,y,alpha=1,beta=1):
    '''
    Potential function for the double well system.
    
    Takes Polarisation and Electric Field
    
    '''
    return beta*x**4 + alpha*x**2 - x*y

def der_F(x,y,alpha=1,beta=1):
    '''
    Derivative of the potential function wrt to x
    
    Takes Polarisation and Electric Field
    
    '''
    return 4*beta*x**3 + 2*alpha*x - y

def der2_F(x,alpha=1,beta=1):
    '''
    Second derivative of the potential function wrt to x
    
    Takes Polarisation and Electric Field
    
    '''
    return 2*alpha + 12*beta*x**2

def findPAPC(e_bias,alpha,beta):
    xminL = sc.optimize.fsolve(der_F,-1,args=(e_bias,alpha,beta))
    xmax = sc.optimize.fsolve(der_F,0,args=(e_bias,alpha,beta))
    xminR = sc.optimize.fsolve(der_F,2,args=(e_bias,alpha,beta))
    return np.array([xminL, xmax, xminR])


def bistable_system(vals,hlev=1,llev=-1,hstate=1,lstate=-1):
    state = lstate
    out_arr = np.zeros_like(vals)
    for i,val in enumerate(vals):
        if state == hstate:
            if val<llev:
                state = lstate
        else:
            if val>hlev:
                state = hstate
        out_arr[i] = state
    return out_arr

# Function for extracting kramers time
def kramerstime_calc(V,P_kramers,high,w,dt,P_thresh=0,switch_counts=[]):
    # Kramers Time Calculation Logic
    mode_counter = bistable_system(V,hlev=-high*0.9,llev=high*0.9)
    mode_count = np.sum(mode_counter)
    v_tracker = np.zeros_like(V,dtype=bool)
    if np.abs(mode_count)/len(mode_counter) <0.5:
        v_tracker[:-1] = np.logical_xor(V[1:]>0,V[:-1]>0)
    else:
        v_tracker[:-1] = np.logical_and(np.logical_not(np.sign(mode_count)*V[1:]>0),np.sign(mode_count)*V[:-1]>0)
    cycles = sum(v_tracker) #total number of cycles
    p_tracker = np.zeros_like(P_kramers,dtype=bool)
    p_tracker[:-1] = np.logical_xor((P_kramers[1:]-P_thresh)>0,(P_kramers[:-1]-P_thresh)>0) #any edge index
    vals =  [int(i) for i in np.where(v_tracker)[0]] + [None] # All the places where voltage switches

    for cycler in range(cycles):
        tmp = np.where(p_tracker[vals[cycler]:vals[cycler+1]])[0] # Selection of when it switches
        if len(tmp) != 0:
            kt_temp = tmp[0]*dt 
            if kt_temp < w*0.9 and kt_temp>0:
                switch_counts.append(kt_temp)
    return switch_counts

def trng_calculator(switch_times_array,bins_array,output_array=[]):
    for thresh_val in bins_array:
        output_array.append(np.sum(switch_times_array<thresh_val)/len(switch_times))
    return output_array

# Functions for fitting
def exponential(x, b):
    return 1- np.exp(-x/b)

# TODO. Pass the variables through the function here
def calculate_kramer_curve(sigma,bias):# Theory Calculation
    R_F = rho*tf/Af # Resistance of the sample
    Dext = sigma**2/(2*R_F*del_F*tf*Af)
    val = findPAPC(bias/tf,alpha,beta)
    P_A = val[0]
    P_C = val[1]
    # print(del_F)
    tau_k = 2*np.pi*rho/(np.sqrt(np.abs(der2_F(P_A,alpha,beta)*der2_F(P_C,alpha,beta))))*np.exp((F(P_C,bias/tf,alpha,beta)-F(P_A,bias/tf,alpha,beta))/Dext)
    return tau_k
#%% Dataset configuration
#----- Script parameters Begin -----#
plt.close('all')
data_gen = 1
histogram_plot = 0
histogram_analysis = 1
read_kramers_csv = 1
new_plot_kramers_time = 1
#----- Script parameters End   -----#
#------ File Parameters to be copied Begin ------ #
write_file_base = '../Datasets/13042025_KramersTime/'
filename = "param.ini"
#------ File Parameters to be copied End ------ #
config = conparse.ConfigParser()
#------ Test Parameter Begin ------ #
config.read(write_file_base+filename)
# Constructing the signal (Copied from the config file)
Ts = float(config["Experiment_Parameters"]["Ts"])# Set a sampling rate
T = float(config["Experiment_Parameters"]["T"]) # total duration
extra_hyst = int(config["Experiment_Parameters"]["Hysteresis Points"]) # total duration
t = np.arange(-extra_hyst*Ts,T+extra_hyst*Ts,Ts)
tr = float(config["Experiment_Parameters"]["tr"])# Reset pulse duration
w = float(config["Experiment_Parameters"]["w"]) # Set pulse duration
high = float(config["Experiment_Parameters"]["high"]) # Reset pulse amplitude
low = np.array(json.loads(config.get("Experiment_Parameters","low")))
T_ = tr + w
offset = float(config["Experiment_Parameters"]["offset"])
sd = np.array(json.loads(config.get("Experiment_Parameters","sd"))) # Set of standard deviations
#------  Test Parameter End ------ #
#------ Physical Parameter Begin------ #
# Physical parameters
tf = float(config['Physical_Parameters']['tf'])
Af = float(config['Physical_Parameters']['Af'])
del_F = float(config['Physical_Parameters']['del_F'])
#------  Physical Parameter End ------ #
# Model Parameters
rho = float(config['Model_Parameter']['rho'])
alpha = float(config['Model_Parameter']['alpha'])
beta  = float(config['Model_Parameter']['beta'])

#%% Preparation of kramers time data 
if data_gen == 1:
    rel_data = []
    sd_i = 2
    times = 2
    for ap in low[:3]:
        switch_times = []
        file_path = write_file_base + f"pulse_{ap}/"+ "data_sd_0_itr0.txt"
        with open(file_path, 'r',encoding='latin-1') as file:
            lines = file.readlines()
            
        V = []
        
        for line in lines[2:]:
            values = line.strip().split('\t')
            V.append(float(values[2]))
        
        V = np.array(V[extra_hyst:-extra_hyst])
        
        for itr in range(times):
            # Extracting a experiment file
            file_path = write_file_base + f"pulse_{ap}/"+ f"data_sd_{sd[sd_i]}_itr{itr}.txt"
            with open(file_path, 'r',encoding='latin-1') as file:
                lines = file.readlines()
            
            t = []
            Vd = []
            P = []
            
            for line in lines[2:]:
                values = line.strip().split('\t')
                t.append(float(values[1]))
                Vd.append(float(values[2]))
                P.append(float(values[3]))
            
            dt = (t[2]-t[1])*1e-3
            # Extracting the portions of the Experiment
            t_kramers = np.array(t[extra_hyst:-extra_hyst])*1e-3
            Vd_kramers = np.array(Vd[extra_hyst:-extra_hyst])
            P_kramers = np.array(P[extra_hyst:-extra_hyst])
            
            # Kramers Time Calculation Logic
            switch_times = kramerstime_calc(V, P_kramers, high, w,dt,switch_counts=switch_times)
        
        switch_times = np.array(switch_times) * 1e6      # µs
        rel_data.append(switch_times)
        
#%% Histogram plot
# Make three plots in separate figures of the histogram
if histogram_plot == 1:
    plt.figure(0)
    plt.clf()
    for i, switch_times in enumerate(rel_data):
        
    
        ax = plt.subplot(1, 1, 1)
        num_bins = 10
    
        # Use logarithmic binning
        bin_edges = np.logspace(np.log10(switch_times.min()),
                                np.log10(switch_times.max()),
                                num_bins)
        
        ax.hist(switch_times, bins=bin_edges, edgecolor='black', alpha=0.7)
        ax.set_xscale('log')
        ax.set_title(f'Histogram for Pulse {low[i]}, SD {sd[sd_i]}')
        ax.set_xlabel('Switch Time (log scale)')
        ax.set_ylabel('Frequency')
        plt.tight_layout()
    plt.savefig("../Figures/histogram.pdf", format="pdf")
#%% TRNG Plot
if histogram_analysis == 1:
    fig4 = plt.figure(4)
    plt.clf()
    fig4.set_size_inches(4, 4, forward=True)
    ax1 = plt.subplot(1, 1, 1)
    
    num_bins = 10
    colors = ['blue', 'green', 'red']  # Add more if rel_data has more items
    labels = ['Bias 1', 'Bias 2', 'Bias 3']
    
    for i, switch_times in enumerate(rel_data):
        bin_edges = np.logspace(np.log10(switch_times.min()*0.2),
                                np.log10(switch_times.max()*5),
                                num_bins+2)
        bin_edges_theory = np.logspace(np.log10(switch_times.min()*0.2),
                                       np.log10(switch_times.max()*5),
                                       num_bins * 100)
    
        trng_data = trng_calculator(switch_times, bin_edges, output_array=[])

        # Plot empirical data
        ax1.scatter(bin_edges, trng_data, label=f'{labels[i]}', color=colors[i], marker='o')
        
        # Fit the exponential model
        try:
            popt, _ = curve_fit(exponential, bin_edges, trng_data, p0=[200], maxfev=5000)
            fitted_b = popt[0]
            fit_curve = exponential(bin_edges_theory, fitted_b)
            ax1.plot(bin_edges_theory, fit_curve,
                     label=f'Theoretical Fit',
                     color=colors[i], linestyle='--')
        except RuntimeError:
            print(f"Fit failed for {labels[i]}")

    ax1.set_xscale('log')
    # ax1.set_title('Normalized Cumulative Distribution of Kramers Time for 0.5V SD')
    ax1.set_xlabel('Time (µs)')
    ax1.grid()
    ax1.set_ylabel('Switching Probability ')
    ax1.legend()
    plt.tight_layout()
    plt.savefig("../Figures/trng_.pdf", format="pdf")

#%%

if read_kramers_csv == 1:
    df_read = pd.read_csv(write_file_base+"kramers_time.csv")
    kramer_time = np.array([df_read[j] for j in ["Bias = " + str(i) + "V" for i in low]],dtype=np.float32)
    sd_true = np.array([df_read[j] for j in ["SD for " + str(i) + "V" for i in low]],dtype=np.float32)


if new_plot_kramers_time == 1:
    fig = plt.figure(2)
    plt.clf()
    fig.set_size_inches(8, 4, forward=True)
    gs = GridSpec(3, 5)
    colors = ['blue', 'green', 'red']
    ax1 = fig.add_subplot(gs[:3,2:])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,1],sharex=ax2)
    ax4 = fig.add_subplot(gs[2,1],sharex=ax2)
    # ax5 = fig.add_subplot(gs[3,:])
    axes_histogram = [ax2,ax3,ax4]
    # annotate_axes(fig)
    # for the kramers time plots
    sigma = np.linspace(0.43,0.7,100)
    Vimprint = -0.3
    for ind in range(len(low[:3])):
        t1 = calculate_kramer_curve(sigma, low[ind]-Vimprint)
        ax1.scatter(sd_true[ind],kramer_time[ind]*1e6,label=f"Experimental Data",color=colors[ind])
        ax1.plot(sigma,t1*1e6,label=f"Theoretical Fit",color=colors[ind])
        
    rect = patches.Rectangle((0.51, 30), 0.03,170, linewidth=1, edgecolor='none', facecolor='orange',alpha=0.2)
    ax1.add_patch(rect)
    ax1.legend()
    ax1.grid()
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.set_xlabel("Standard Deviation(V)")
    ax1.set_ylabel("Time($\mu s$)")
    
    
    
    
    # for the histogram plots
    for i, switch_times in enumerate(rel_data):
        num_bins = 15
    
        # Use logarithmic binning
        bin_edges = np.logspace(np.log10(switch_times.min()),
                                np.log10(switch_times.max()),
                                num_bins)
    
        axes_histogram[i].hist(switch_times, bins=bin_edges, edgecolor='black', alpha=0.7,color = colors[i])
        axes_histogram[i].set_xscale('log')
        #xes_histogram[i].plot(bin_edges,40e12*bin_edges*np.exp(-bin_edges/120e-6))
        # axes_histogram[i].set_title(f'Histogram for Pulse {low[i]}, SD {sd[sd_i]}')
        # axes_histogram[i].set_xlabel('Switch Time (log scale)')
        # axes_histogram[i].set_ylabel('Frequency')
        if True:
            axes_histogram[i].tick_params(axis='y', labelbottom=False) 
            # print(i)
    
    plt.tight_layout()
    plt.show()
    plt.savefig("../Figures/kramers_time.pdf", format="pdf")
    

    
    