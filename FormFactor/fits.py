#!/usr/bin/python3
import numpy as np
import math
import cmath
import random
from scipy import optimize
import matplotlib.pyplot as plt


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
            for config in range(n_configs):
                cov_mat[x,y] += (corr[config,x] - av_corr[x])*(corr[config,y] - av_corr[y])
            cov_mat[x,y] /= float(n_configs)*float(n_configs-1.)
    w_mat = np.linalg.inv(cov_mat)
    return w_mat


def chi2(param, av_corr, w_mat, fit_func, Nt, nmin, nmax):
    x2 = 0.
    for x in range(nmin, nmax+1):
        for y in range(nmin, nmax+1):
            x2 += (av_corr[x] - fit_func(x, Nt, param))*w_mat[x,y]*(av_corr[y] - fit_func(y, Nt, param))
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
    nmin = nmin
    nmax = nmax
    res = optimize.minimize(chi2, par_guess, method='Nelder-Mead', args=(av_corr, w_mat, fit_func, Nt, nmin, nmax))
    if(res.success):
        print('\033[1;32mMinimization succeeded\033[0m')
        print('A0 = '+str(res.x[0])+', m = '+str(res.x[1]))
        return res.x
    else:
        print('\033[1;31mMinimization failed\033[0m')
        print(res.message)
        return [0,0]
