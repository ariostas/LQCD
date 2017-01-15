#!/usr/bin/python3
import numpy as np
import math
import cmath
import random
import matplotlib.pyplot as plt
import fits
import multiprocessing
import lqcdfunc

# file_prefix = '/home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450/'
file_prefix = '/home/arios/Documents/LQCDConfigs/wil_16_64_aniso/'

sources = ['DG0_1', 'DG1_1', 'DG1_1', 'DG2_1', 'DG2_1']
sinks = ['DG0_1', 'DG1_1', 'DG1_1', 'DG2_1', 'DG2_1']

size_s = 16 # spatial size
size_t = 64 # temporal size
t0 = 0 # time used for generalized eigenvalue problem # t = 0 for both proton and rho
t1 = 18 # time at which variational source is picked, and used for S/N optimization # t = 18 for proton and t = 10 for rho
t_sink = 16 # sink time (for 3ptfns)

proton_3ptfn_t_corr = lqcdfunc.read_file(data_type='3ptfn', seqsource='NUCL_D_UNPOL', g=7, sources=sources, sinks=sinks, q=(0,0,0), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
proton_2ptfn_srcp_corr = lqcdfunc.read_file(data_type='2ptfn', had='proton', pf=(0,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)
proton_2ptfn_snkp_corr = lqcdfunc.read_file(data_type='2ptfn', had='proton', pf=(0,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)

proton_3ptfn_t_mat = lqcdfunc.construct_matrices(proton_3ptfn_t_corr)
proton_2ptfn_srcp_mat = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr)
proton_2ptfn_snkp_mat = lqcdfunc.construct_matrices(proton_2ptfn_snkp_corr)

proton_3ptfn_t_mat = lqcdfunc.make_hermitian(proton_3ptfn_t_mat)
proton_2ptfn_srcp_mat = lqcdfunc.make_hermitian(proton_2ptfn_srcp_mat)
proton_2ptfn_snkp_mat = lqcdfunc.make_hermitian(proton_2ptfn_snkp_mat)

proton_ff_t, proton_fferr_t = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, t_sink, t0, t1)
proton_mass, proton_masserr = lqcdfunc.find_gnd_masses(proton_2ptfn_snkp_mat, t0, t1)


xlab = np.arange(0, size_t/2, 1)
type_labels = ['Point', '$G1$', '$\\nabla^2 G1$', '$G2$', '$\\nabla^2 G2$', 'Var', 'Var + S/N', 'S/N']
plt.figure(1, figsize=(16, 12))
for x in range(5):
    plt.subplot(511+x)
    plt.errorbar(xlab, np.resize(proton_mass[x], size_t/2), yerr=np.resize(proton_masserr[x], size_t/2), color='b', ecolor='b', fmt='o', capsize=2)
    plt.ylabel('$m_{eff}$', fontsize=16)
    plt.xlabel('$ \Delta t $ (%s)'  % type_labels[x], fontsize=16)
    plt.xlim(-0.5, size_t/2+0.5)
    # plt.yscale('log', nonposy='clip')
    # plt.ylim(0, 2)
    # xp = np.linspace(fitrange[x][0], fitrange[x][1], 100)
    # yp = np.zeros(len(xp))
    # for i in range(len(xp)):
    #     yp[i] = masses[x]
    # plt.plot(xp, yp, color='r')

type_labels = ['Var', 'Var+S/N', 'S/N']
colors = ['r', 'b', 'g']
markers = ['^', 's', 'o']
plt.figure(2, figsize=(10, 7))
for x in range(3):
    ax = plt.subplot(311+x)
    plt.errorbar(xlab, np.resize(proton_mass[x+5], size_t/2), yerr=np.resize(proton_masserr[x+5], size_t/2), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2)
    plt.ylabel('$m_{eff}$', fontsize=20)
    # plt.xlabel('$ \Delta t $ (%s)'  % type_labels[x], fontsize=16)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, size_t/2+0.5)
    plt.ylim(0, 1.4)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 0.2))
    if(x == 2):
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plt.xlabel('$ \Delta t $', fontsize=20)
    else:
        ax.set_xticklabels( () )
    if(x == 0):
        plt.suptitle('Proton', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)


plt.figure(3, figsize=(10, 7))
for x in range(3):
    ax = plt.subplot(311+x)
    plt.errorbar(np.resize(xlab, t_sink), np.resize(proton_ff_t[x+5], t_sink), yerr=np.resize(proton_fferr_t[x+5], t_sink), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2)
    plt.ylabel('$R$', fontsize=20)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, t_sink+0.5)
    plt.ylim(0.0, 0.3)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 0.1))
    if(x == 2):
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plt.xlabel('$ \\tau $', fontsize=20)
    else:
        ax.set_xticklabels( () )
    if(x == 0):
        plt.suptitle('Proton ($\gamma_{4}\gamma_{5}$ charge)', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

type_labels = ['Point', '$G1$', '$\\nabla^2 G1$', '$G2$', '$\\nabla^2 G2$', 'Var', 'Var + S/N', 'S/N']
colors = ['c', 'b', 'g', 'k', 'm', 'r']
markers = ['^', 's', 'o', '8', 'p', '*']
plt.figure(4, figsize=(10, 10))
for x in range(5):
    ax = plt.subplot(511+x)
    plt.errorbar(np.resize(xlab, t_sink), np.resize(proton_ff_t[x], t_sink), yerr=np.resize(proton_fferr_t[x], t_sink), color='b', ecolor='b', fmt='o', capsize=2)
    plt.ylabel('$R$', fontsize=20)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, t_sink+0.5)
    # plt.ylim(0.0, 0.25)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 0.1))
    if(x == 4):
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plt.xlabel('$ \\tau $', fontsize=20)
    else:
        ax.set_xticklabels( () )
    if(x == 0):
        plt.suptitle('Proton axial charge', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

plt.show()


