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
file_prefix = '/home/arios/Documents/LQCDConfigs/wil_16_60_aniso_cluster/rho/'

sources = ['DG0_1', 'DG1_1', 'DG1_1', 'DG2_1', 'DG2_1']
sinks = ['DG0_1', 'DG1_1', 'DG1_1', 'DG2_1', 'DG2_1']

size_s = 16 # spatial size
size_t = 60 # temporal size
t0 = 0 # time used for generalized eigenvalue problem # t = 0 for both proton and rho
t1 = 10 # time at which variational source is picked, and used for S/N optimization # t = 18 for proton and t = 10 for rho
t_sink = 21 # sink time (for 3ptfns)
n_boot = 100

proton_3ptfn_t_corr1 = lqcdfunc.read_file(data_type='3ptfn', seqsource='a0-rho_x_1', i=1, g=8, sources=sources, sinks=sinks, q=(1,0,0), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
proton_3ptfn_t_corr2 = lqcdfunc.read_file(data_type='3ptfn', seqsource='a0-rho_x_1', i=1, g=8, sources=sources, sinks=sinks, q=(-1,0,0), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
# proton_3ptfn_t_corr3 = lqcdfunc.read_file(data_type='3ptfn', seqsource='a0-rho_x_1', i=1, g=8, sources=sources, sinks=sinks, q=(0,1,0), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
# proton_3ptfn_t_corr4 = lqcdfunc.read_file(data_type='3ptfn', seqsource='a0-rho_x_1', i=1, g=8, sources=sources, sinks=sinks, q=(0,-1,0), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
# proton_3ptfn_t_corr5 = lqcdfunc.read_file(data_type='3ptfn', seqsource='a0-rho_x_1', i=1, g=8, sources=sources, sinks=sinks, q=(0,0,1), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
# proton_3ptfn_t_corr6 = lqcdfunc.read_file(data_type='3ptfn', seqsource='a0-rho_x_1', i=1, g=8, sources=sources, sinks=sinks, q=(0,0,-1), pf=(0,0,0), current='nonlocal', file_prefix=file_prefix)
proton_2ptfn_srcp_corr1 = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(-1,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)
proton_2ptfn_srcp_corr2 = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(1,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)
# proton_2ptfn_srcp_corr3 = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(0,-1,0), sources=sources, sinks=sinks, file_prefix=file_prefix)
# proton_2ptfn_srcp_corr4 = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(0,1,0), sources=sources, sinks=sinks, file_prefix=file_prefix)
# proton_2ptfn_srcp_corr5 = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(0,0,-1), sources=sources, sinks=sinks, file_prefix=file_prefix)
# proton_2ptfn_srcp_corr6 = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(0,0,1), sources=sources, sinks=sinks, file_prefix=file_prefix)
proton_2ptfn_snkp_corr = lqcdfunc.read_file(data_type='2ptfn', had='rho_x', pf=(0,0,0), sources=sources, sinks=sinks, file_prefix=file_prefix)

proton_3ptfn_t_mat1 = lqcdfunc.construct_matrices(proton_3ptfn_t_corr1)
proton_3ptfn_t_mat2 = lqcdfunc.construct_matrices(proton_3ptfn_t_corr2)
# proton_3ptfn_t_mat3 = lqcdfunc.construct_matrices(proton_3ptfn_t_corr3)
# proton_3ptfn_t_mat4 = lqcdfunc.construct_matrices(proton_3ptfn_t_corr4)
# proton_3ptfn_t_mat5 = lqcdfunc.construct_matrices(proton_3ptfn_t_corr5)
# proton_3ptfn_t_mat6 = lqcdfunc.construct_matrices(proton_3ptfn_t_corr6)
proton_2ptfn_srcp_mat1 = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr1)
proton_2ptfn_srcp_mat2 = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr2)
# proton_2ptfn_srcp_mat3 = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr3)
# proton_2ptfn_srcp_mat4 = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr4)
# proton_2ptfn_srcp_mat5 = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr5)
# proton_2ptfn_srcp_mat6 = lqcdfunc.construct_matrices(proton_2ptfn_srcp_corr6)
proton_2ptfn_snkp_mat = lqcdfunc.construct_matrices(proton_2ptfn_snkp_corr)

# proton_3ptfn_t_mat = lqcdfunc.make_hermitian(proton_3ptfn_t_mat)
# proton_2ptfn_srcp_mat = lqcdfunc.make_hermitian(proton_2ptfn_srcp_mat)
# proton_2ptfn_snkp_mat = lqcdfunc.make_hermitian(proton_2ptfn_snkp_mat)

proton_ff_t, proton_fferr_t = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat1, proton_2ptfn_srcp_mat1, proton_2ptfn_snkp_mat, t_sink, t0, t1, n_boot=n_boot)
proton_ff_t2, proton_fferr_t2 = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat2, proton_2ptfn_srcp_mat2, proton_2ptfn_snkp_mat, t_sink, t0, t1, n_boot=n_boot)
# proton_ff_t3, proton_fferr_t3 = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat3, proton_2ptfn_srcp_mat3, proton_2ptfn_snkp_mat, t_sink, t0, t1, n_boot=n_boot)
# proton_ff_t4, proton_fferr_t4 = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat4, proton_2ptfn_srcp_mat4, proton_2ptfn_snkp_mat, t_sink, t0, t1, n_boot=n_boot)
# proton_ff_t5, proton_fferr_t5 = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat5, proton_2ptfn_srcp_mat5, proton_2ptfn_snkp_mat, t_sink, t0, t1, n_boot=n_boot)
# proton_ff_t6, proton_fferr_t6 = lqcdfunc.find_ff_with2ptfn(proton_3ptfn_t_mat6, proton_2ptfn_srcp_mat6, proton_2ptfn_snkp_mat, t_sink, t0, t1, n_boot=n_boot)
proton_mass, proton_masserr = lqcdfunc.find_gnd_masses(proton_2ptfn_snkp_mat, t0, t1, n_boot=n_boot)


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
        plt.suptitle('Rho', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)


plt.figure(3, figsize=(10, 7))
for x in range(3):
    ax = plt.subplot(311+x)
    plt.errorbar(np.resize(xlab, t_sink), np.resize(proton_ff_t[x+5], t_sink), yerr=np.resize(proton_fferr_t[x+5], t_sink), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2)
    plt.ylabel('$R$', fontsize=20)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, t_sink+0.5)
    plt.ylim(0.3, 0.9)
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
        plt.suptitle('Rho ($\gamma_4$-current)', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

type_labels = ['Point', '$G1$', '$\\nabla^2 G1$', '$G2$', '$\\nabla^2 G2$', 'Var', 'Var + S/N', 'S/N']
colors = ['c', 'b', 'g', 'k', 'm', 'r']
markers = ['^', 's', 'o', '8', 'p', '*']
# plt.figure(4, figsize=(16, 12))
# for x in range(5):
#     plt.subplot(511+x)
#     plt.errorbar(np.resize(xlab, t_sink), np.resize(proton_ff_t[x], t_sink), yerr=np.resize(proton_fferr_t[x], t_sink), color='b', ecolor='b', fmt='o', capsize=2)
#     plt.ylabel('$R$', fontsize=16)
#     plt.xlabel('$ \\tau $ (%s)'  % type_labels[x], fontsize=16)
#     plt.xlim(-0.5, t_sink+0.5)

plt.figure(4, figsize=(10, 6))
for x in range(5):
    plt.errorbar(np.resize(xlab, t_sink)+x*0.1, np.resize(proton_ff_t[x], t_sink), yerr=np.resize(proton_fferr_t[x], t_sink), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2, label=type_labels[x])
plt.errorbar(np.resize(xlab, t_sink)+0.5, np.resize(proton_ff_t[5], t_sink), yerr=np.resize(proton_fferr_t[5], t_sink), color=colors[5], ecolor=colors[5], fmt=markers[5], capsize=2, label='Variational')
ax = plt.subplot(111)
plt.ylabel('$R$', fontsize=20)
plt.xlabel('$ \\tau $', fontsize=20)
plt.xlim(-0.5, t_sink+0.5)
plt.ylim(0.3, 1.1)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode='expand', borderaxespad=0.)
plt.suptitle('Variational vs diagonal correlators', fontsize=30)
plt.text(0.45, 0.82, 'Rho $\gamma_4$-current q=(1,0,0)', transform=ax.transAxes, fontsize=20)
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
start, end = ax.get_ylim()
ax.yaxis.set_ticks(np.arange(start, end, 0.1))
plt.subplots_adjust(left=0.10, right=0.96, top=0.76, bottom=0.11, wspace=0.22, hspace=0.12)

plt.figure(5, figsize=(10, 6))
for x in range(5):
    plt.errorbar(np.resize(xlab, t_sink)+x*0.1, np.resize(proton_ff_t[x], t_sink), yerr=np.resize(proton_fferr_t[x], t_sink), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2, label=type_labels[x])
plt.errorbar(np.resize(xlab, t_sink)+0.5, np.resize(proton_ff_t[7], t_sink), yerr=np.resize(proton_fferr_t[7], t_sink), color=colors[5], ecolor=colors[5], fmt=markers[5], capsize=2, label='S/N')
ax = plt.subplot(111)
plt.ylabel('$R$', fontsize=20)
plt.xlabel('$ \\tau $', fontsize=20)
plt.xlim(-0.5, t_sink+0.5)
plt.ylim(0.3, 1.1)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode='expand', borderaxespad=0.)
plt.suptitle('S/N optimized vs diagonal correlators', fontsize=30)
plt.text(0.45, 0.82, 'Rho $\gamma_4$-current q=(1,0,0)', transform=ax.transAxes, fontsize=20)
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
start, end = ax.get_ylim()
ax.yaxis.set_ticks(np.arange(start, end, 0.1))
plt.subplots_adjust(left=0.10, right=0.96, top=0.76, bottom=0.11, wspace=0.22, hspace=0.12)

plt.figure(5, figsize=(10, 10))
for x in range(5):
    ax = plt.subplot(511+x)
    plt.errorbar(np.resize(xlab, t_sink), np.resize(proton_ff_t[x], t_sink), yerr=np.resize(proton_fferr_t[x], t_sink), color='b', ecolor='b', fmt='o', capsize=2)
    plt.ylabel('$R$', fontsize=20)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, t_sink+0.5)
    plt.ylim(0.3, 0.9)
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
        plt.suptitle('Rho ($\gamma_4$-current)', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

plt.figure(6, figsize=(10, 6))
for x in range(5):
    plt.errorbar(np.resize(xlab, size_t/2)+x*0.1, np.resize(proton_mass[x], size_t/2), yerr=np.resize(proton_masserr[x], size_t/2), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2, label=type_labels[x])
plt.errorbar(np.resize(xlab, size_t/2)+0.5, np.resize(proton_mass[7], size_t/2), yerr=np.resize(proton_masserr[7], size_t/2), color=colors[5], ecolor=colors[5], fmt=markers[5], capsize=2, label='S/N')
ax = plt.subplot(111)
plt.ylabel('$m_{eff}$', fontsize=20)
plt.xlabel('$ \Delta t $', fontsize=20)
plt.xlim(-0.5, size_t/2+0.5)
# plt.ylim(0.3, 1.1)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode='expand', borderaxespad=0.)
plt.suptitle('S/N optimized vs diagonal correlators', fontsize=30)
plt.text(0.45, 0.82, 'Rho 2ptfn', transform=ax.transAxes, fontsize=20)
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
start, end = ax.get_ylim()
ax.yaxis.set_ticks(np.arange(start, end, 0.2))
plt.subplots_adjust(left=0.10, right=0.96, top=0.76, bottom=0.11, wspace=0.22, hspace=0.12)

plt.figure(7, figsize=(10, 6))
for x in range(5):
    plt.errorbar(np.resize(xlab, size_t/2)+x*0.1, np.resize(proton_mass[x], size_t/2), yerr=np.resize(proton_masserr[x], size_t/2), color=colors[x], ecolor=colors[x], fmt=markers[x], capsize=2, label=type_labels[x])
plt.errorbar(np.resize(xlab, size_t/2)+0.5, np.resize(proton_mass[5], size_t/2), yerr=np.resize(proton_masserr[5], size_t/2), color=colors[5], ecolor=colors[5], fmt=markers[5], capsize=2, label='Variational')
ax = plt.subplot(111)
plt.ylabel('$m_{eff}$', fontsize=20)
plt.xlabel('$ \Delta t $', fontsize=20)
plt.xlim(-0.5, size_t/2+0.5)
# plt.ylim(0.3, 1.1)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode='expand', borderaxespad=0.)
plt.suptitle('Variational vs diagonal correlators', fontsize=30)
plt.text(0.45, 0.82, 'Rho 2ptfn', transform=ax.transAxes, fontsize=20)
for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
start, end = ax.get_ylim()
ax.yaxis.set_ticks(np.arange(start, end, 0.2))
plt.subplots_adjust(left=0.10, right=0.96, top=0.76, bottom=0.11, wspace=0.22, hspace=0.12)

plt.figure(8, figsize=(10, 10))
for x in range(5):
    ax = plt.subplot(511+x)
    plt.errorbar(np.resize(xlab, size_t/2), np.resize(proton_mass[x], size_t/2), yerr=np.resize(proton_masserr[x], size_t/2), color='b', ecolor='b', fmt='o', capsize=2)
    plt.ylabel('$m_{eff}$', fontsize=20)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, size_t/2+0.5)
    plt.ylim(0, 1.0)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 0.2))
    if(x == 4):
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plt.xlabel('$ \Delta t $', fontsize=20)
    else:
        ax.set_xticklabels( () )
    if(x == 0):
        plt.suptitle('Rho', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

type_labels = ['q=(1,0,0)', 'q=(-1,0,0)', 'q=(0,1,0)', 'q=(0,-1,0)', 'q=(0,0,1)', 'q=(0,0,-1)']
colors = ['r', 'b', 'g']
markers = ['^', 's', 'o']
dat = [proton_ff_t, proton_ff_t2]#, proton_ff_t3, proton_ff_t4, proton_ff_t5, proton_ff_t6]
err = [proton_fferr_t, proton_fferr_t2]#, proton_fferr_t3, proton_fferr_t4, proton_fferr_t5, proton_fferr_t6]
plt.figure(9, figsize=(16, 10))
for x in range(2):
    ax = plt.subplot(211+x)
    plt.errorbar(np.resize(xlab, t_sink), np.resize(dat[x][0], t_sink), yerr=np.resize(err[x][0], t_sink), color='b', ecolor='b', fmt='o', capsize=2)
    plt.ylabel('$R$', fontsize=20)
    # plt.xlabel('$ \Delta t $ (%s)'  % type_labels[x], fontsize=16)
    plt.text(0.45, 0.82, type_labels[x], transform=ax.transAxes, fontsize=20)
    plt.xlim(-0.5, t_sink+0.5)
    plt.ylim(0.3, 1.0)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 0.1))
    if(x==1):#if(x == 4 or x == 5):
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plt.xlabel('$ \\tau $', fontsize=20)
    else:
        ax.set_xticklabels( () )
    if(x == 0):
        plt.suptitle('Rho ($\gamma_4$-current)', fontsize=30)
plt.subplots_adjust(left=0.08, right=0.97, top=0.92, bottom=0.10, wspace=0.22, hspace=0.12)

plt.show()
