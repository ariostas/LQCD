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

# t_0 = 0 and t_1 = 6 work good
t0 = 0 # time used for generalized eigenvalue problem
t1 = 8 # time at which variational source is picked, and used for S/N optimization
t_sink = 16

proton_3ptfn_t_corr = lqcdfunc.read_file(data_type='3ptfn', g=8, sources=sources, sinks=sinks, q=(1,0,0), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
proton_2ptfn_srcp_corr = lqcdfunc.read_file(data_type='2ptfn', had='proton', pf=(-1,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)
proton_2ptfn_snkp_corr = lqcdfunc.read_file(data_type='2ptfn', had='proton', pf=(0,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)

proton_3ptfn_t_mat = lqcdfunc.construct_matrices(proton_3ptfn_t_corr)
proton_2ptfn_srcp_mat = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr)
proton_2ptfn_snkp_mat = lqcdfunc.construct_matrices(proton_2ptfn_snkp_corr)
proton_rec3ptfn_mat = lqcdfunc.construct_rec_matrices(proton_3ptfn_t_corr, 8, 8, 16)

proton_3ptfn_t_mat = lqcdfunc.make_hermitian(proton_3ptfn_t_mat)
proton_2ptfn_srcp_mat = lqcdfunc.make_hermitian(proton_2ptfn_srcp_mat)
proton_2ptfn_snkp_mat = lqcdfunc.make_hermitian(proton_2ptfn_snkp_mat)

proton_ff_t, proton_fferr_t = lqcdfunc.find_rec_ff(proton_3ptfn_t_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, proton_rec3ptfn_mat, t_sink, t0, t1)

xlab = np.arange(0, t_sink, 1)
type_labels = ['Var', 'Var+S/N', 'S/N']
colors = ['r', 'b', 'g']
markers = ['^', 's', 'o']
plt.figure(3, figsize=(10, 7))
for x in range(3):
    ax = plt.subplot(311+x)
    plt.errorbar(np.resize(xlab, t_sink), np.resize(proton_ff_t[x+5], t_sink), yerr=np.resize(proton_fferr_t[x+5], t_sink), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2)
    plt.ylabel('$R$', fontsize=20)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, t_sink+0.5)
    # plt.ylim(0.5, 1.0)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    # start, end = ax.get_ylim()
    # ax.yaxis.set_ticks(np.arange(start, end, 0.1))
    if(x == 2):
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plt.xlabel('$ \\tau $', fontsize=20)
    else:
        ax.set_xticklabels( () )
    if(x == 0):
        plt.suptitle('Proton ($\gamma_5$-current)', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

plt.show()
