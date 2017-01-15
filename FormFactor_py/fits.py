#!/usr/bin/python3
import numpy as np
import math
import cmath
import random
from scipy import optimize
import matplotlib.pyplot as plt
import multiprocessing


def average(corr):
    n_configs = len(corr)
    av = corr[0]
    for x in range(1, n_configs):
        av += corr[x]
    av /= float(n_configs)
    return av


def exp_func(nt, Nt, par):
    return par[0]*math.exp(-par[1]*nt)


def cosh_func(nt, Nt, par):
    return par[0]*math.cosh((nt-Nt/2.)*par[1])


def sinh_func(nt, Nt, par):
    return par[0]*math.sinh((nt-Nt/2.)*par[1])


def weight_mat(corr):
    n_configs, Nt = corr.shape
    av_corr = average(corr)
    cov_mat = np.matrix(np.zeros((Nt,Nt)))
    for x in range(Nt):
        for y in range(Nt):
            # if(x != y):
            #     continue
            for config in range(n_configs):
                cov_mat[x,y] += (corr[config,x] - av_corr[x])*(corr[config,y] - av_corr[y])
            cov_mat[x,y] /= float(n_configs)*float(n_configs-1.)
    w_mat = np.linalg.inv(cov_mat)
    # print(cov_mat)
    # print(w_mat)
    return w_mat


def chi2(param, av_corr, w_mat, fit_func, Nt, nmin, nmax):
    x2 = 0.
    for x in range(nmin, nmax+1):
        for y in range(nmin, nmax+1):
            x2 += (av_corr[x] - fit_func(x, Nt, param))*w_mat[x,y]*(av_corr[y] - fit_func(y, Nt, param))
    # x2 = x2/float(nmax - nmin - 1)
    return x2


def fit_corr(corr, nmin, nmax, par_guess, func='exp'):
    corr = np.array(corr).real
    if(func == 'exp'):
        fit_func = exp_func
    elif(func == 'cosh'):
        fit_func = cosh_func
    elif(func == 'sinh'):
        fit_func = sinh_func
    else:
        print('Function type not defined, setting it to exp')
        fit_func = exp_func
    av_corr = average(corr)
    Nt = len(av_corr)
    w_mat = weight_mat(corr)
    res = optimize.minimize(chi2, par_guess, method='Nelder-Mead', args=(av_corr, w_mat, fit_func, Nt, nmin, nmax))
    if(res.success):
        # print('\033[1;32mMinimization succeeded\033[0m')
        # print('A0 = '+str(res.x[0])+', m = '+str(res.x[1]))
        chi2_min = chi2(res.x, av_corr, w_mat, fit_func, Nt, nmin, nmax)
        # print('min chi2 = '+str(chi2_min))
        return res.x
    else:
        print('\033[1;31mMinimization failed\033[0m')
        print(res.message)
        print('A0 = '+str(res.x[0])+', m = '+str(res.x[1]))
        chi2_min = chi2(res.x, av_corr, w_mat, fit_func, Nt, nmin, nmax)
        print('min chi2 = '+str(chi2_min))
        return [0,0]


def fit_correrr(corr, nmin, nmax, par_guess, func='exp'):
    corr = np.array(corr).real
    if(func == 'exp'):
        fit_func = exp_func
    elif(func == 'cosh'):
        fit_func = cosh_func
    elif(func == 'sinh'):
        fit_func = sinh_func
    else:
        print('Function type not defined, setting it to exp')
        fit_func = exp_func
    av_corr = average(corr)
    Nt = len(av_corr)
    w_mat = weight_mat(corr)
    res = optimize.minimize(chi2, par_guess, method='Nelder-Mead', args=(av_corr, w_mat, fit_func, Nt, nmin, nmax))
    if(res.success):
        # print('\033[1;32mMinimization succeeded\033[0m')
        # print('A0 = '+str(res.x[0])+', m = '+str(res.x[1]))
        chi2_min = chi2(res.x, av_corr, w_mat, fit_func, Nt, nmin, nmax)
        # print('min chi2 = '+str(chi2_min))
        step = 0.000005
        m_low = res.x[1]
        m_high = res.x[1]
        amp = res.x[0]
        while(chi2((amp,m_low), av_corr, w_mat, fit_func, Nt, nmin, nmax) - chi2_min < 1):
            m_low = m_low - step
            # print(str(chi2(param, av_corr, w_mat, fit_func, Nt, nmin, nmax))+'  '+str(m_low)+'  '+str(res.x[1]))
        while(chi2((amp,m_high), av_corr, w_mat, fit_func, Nt, nmin, nmax) - chi2_min < 1):
            m_high = m_high + step
        # print(str(m_high)+'  '+str(m_low))
        return res.x, (m_high - m_low)/2., chi2_min/float(nmax - nmin - 1)
    else:
        print('\033[1;31mMinimization failed\033[0m')
        print(res.message)
        return [0,0], 0, 0

def scan_fit(corr, par_guess, res, name, nmin=10, nmax=17, func='exp'):
    print('Scanning fit range...')
    n_configs, t_size = corr.shape
    n_low = []
    n_high = []
    chi_low = []
    chi_high = []
    for x in range(nmax-1):
        temp1, temp2, temp3 = fit_correrr(corr, nmax-2-x, nmax, par_guess, func)
        n_low.append(nmax-2-x)
        chi_low.append(temp3)
        if(temp3 > 10):
            break
    for x in range(t_size-nmin-2):
        temp1, temp2, temp3 = fit_correrr(corr, nmin, nmin+2+x, par_guess, func)
        n_high.append(nmin+2+x)
        chi_high.append(temp3)
        if(temp3 > 10):
            break
    # print(chi_low)
    # print(chi_high)
    n_min = 0
    n_max = 0
    for x in range(len(n_low)):
        if(chi_low[x] < 6):
            n_min = n_low[x]
    for x in range(len(n_high)):
        if(chi_high[x] < 3 and n_high[x] < 32):
            n_max = n_high[x]
    all_mass = []
    all_masserr = []
    # print(str(n_min)+'    '+str(n_max))
    for x in range(3):
        for y in range(3):
            temp1, temp2, temp3 = fit_correrr(corr, n_min-1+x, n_max-1+y, par_guess, func)
            all_mass.append(temp1[1])
            all_masserr.append(temp2)
    av_mass = average(all_mass)
    av_masserr = average(all_masserr)
    std_mass = 0
    if(len(all_mass) != 9):
        print('Error: Wrong length')
    for x in range(len(all_mass)):
        std_mass += (all_mass[x]-av_mass)**2
    std_mass = 1./9.*std_mass**(1./2.)
    # print('Mass = '+str(av_mass)+' +- '+str(av_masserr)+' +- '+str(std_mass))
    res.put([name, av_mass, av_masserr, std_mass, n_min, n_max])
    # return av_mass, av_masserr+std_mass, n_min, n_max
