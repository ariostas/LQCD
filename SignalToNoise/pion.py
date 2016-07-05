import numpy as np
import scipy.optimize as opt
import math
import cmath
import random
import matplotlib.pyplot as plt

random.seed()
FILE_NAMES = ["" for x in range(25)]
for x in range(1, 6):
    for y in range(1, 6):
        FILE_NAMES[(x-1)*5+(y-1)] = '../SignalToNoise/PaperData/proton_SrcDG' + str(x) + '_SnkDG' + str(y) + '_Interp4.dat'


def read_file(names):
    f_counter = 0
    n_files = len(names)
    data = np.zeros((1, n_files, 1), dtype=np.complex128)
    for file in names:
        if(f_counter < n_files-1):
            print('Reading correlator data... ' + str(math.floor(f_counter/n_files*100)) + '%', end='\r')
        else:
            print('Reading correlator data... \033[1;32msuccess\033[0m')
        f = open(file, 'r')
        counter = -1
        n_configs = 0
        t_size = 0
        for line in f:
            line = line.strip()
            col = line.split()
            if(counter == -1):
                n_configs = int(col[0])
                t_size = int(col[1])
                if(f_counter == 0):
                    data = np.zeros((int(col[0]), n_files , int(col[1])), dtype=np.complex128)
            else:
                data[math.floor(counter/t_size), f_counter, counter%t_size] = float(col[1]) + float(col[2])*1.0j
            counter += 1
        if(counter != n_configs*t_size):
            print('Error: File seems to be corrupted')
        f_counter += 1
    return data


def average(corr):
    n_configs, n_types, t_size = corr.shape
    av = np.zeros((n_types, t_size), dtype=np.complex128)
    for x in range(0, n_configs):
        av += corr[x]
    av /= float(n_configs)
    return av


def eff_masses(corr):
    t_size = corr.size
    masses = np.zeros(t_size)
    for x in range(0, t_size):
        r = corr[x].real/corr[(x+1)%t_size].real
        m = 0
        if(r > 0.):
            m = math.log(r)
        if(x > t_size/2):
            m = -m
        masses[x] = m
    return masses


def make_matrices(corr):
    n_types, t_size = corr.shape
    mat_size = math.sqrt(n_types)
    if(mat_size%1 != 0):
        print("Matrix size error")
    matrices = [np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128)) for i in range(t_size)]
    for x in range(t_size):
        temp_arr = np.zeros((mat_size, mat_size), dtype=np.complex128)
        for y in range(n_types):
            temp_arr[math.floor(y/mat_size), y%mat_size] = corr[y, x]
        matrices[x] = np.matrix(temp_arr)
    return matrices


def find_eigsys(mat):
    t_size = len(mat)
    mat_size, temp = mat[0].shape
    evals = np.zeros((mat_size, t_size), dtype=np.complex128)
    evecs = [[np.matrix(np.zeros((mat_size, 1)), dtype=np.complex128) for n in range(t_size)] for m in range(mat_size)]
    init = mat[0]
    init2 = mat[3]
    init = np.linalg.inv(init)
    init2 = np.linalg.inv(init2)
    for x in range(t_size):
        temp_evals, temp_evecs = np.linalg.eig(init * mat[x])
        temp_evals2, temp_evecs2 = np.linalg.eig(init2 * mat[x])
        order = np.argsort(temp_evals)
        order2 = np.argsort(temp_evals2)
        for y in range(mat_size):
            evals[y, x] = temp_evals[order[mat_size-1-y]]
            evecs[y][x] = temp_evecs2[:,order2[mat_size-1-y]]
    return evals, evecs


def find_sink(corr, source, mat, ts):
    mat_size, temp = mat[0].shape
    n_configs, n_types, t_size = corr.shape
    sink = np.matrix(np.zeros((mat_size, 1), dtype=np.complex128))
    sigma2 = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
    for x in range(n_configs):
        c = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
        for y in range(n_types):
            c[math.floor(y/mat_size), y%mat_size] = corr[x, y, ts]
        sigma2 += c * source * source.H * c.H
    sigma2 /= n_configs
    inv_sigma2 = np.linalg.inv(sigma2)
    a = (source.H * mat[ts].H * inv_sigma2 * inv_sigma2 * mat[ts] * source)[0,0]
    a = 1.0/cmath.sqrt(a)
    sink = a * inv_sigma2 * mat[ts] * source
    sink /= cmath.sqrt((sink.H * sink)[0,0])
    return sink


def find_opt_sn(source_g, sink_g, corr, mat, ts):
    mat_size, temp = mat[0].shape
    n_configs, n_types, t_size = corr.shape
    n_iter = 20
    source = source_g
    sink = sink_g
    for n in range(n_iter):
        sigma2_source = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
        sigma2_sink = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
        for x in range(n_configs):
            c = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
            for y in range(n_types):
                c[math.floor(y/mat_size), y%mat_size] = corr[x, y, ts]
            sigma2_source += c * source * source.H * c.H
            sigma2_sink += c.H * sink * sink.H * c
        sigma2_source /= n_configs
        sigma2_sink /= n_configs
        inv_sigma2_source = np.linalg.inv(sigma2_source)
        inv_sigma2_sink = np.linalg.inv(sigma2_sink)
        a_source = (source.H * mat[ts].H * inv_sigma2_source * inv_sigma2_source * mat[ts] * source)[0,0]
        a_sink = (sink.H * mat[ts] * inv_sigma2_sink * inv_sigma2_sink * mat[ts].H * sink)[0,0]
        a_source = 1.0/cmath.sqrt(a_source)
        a_sink = 1.0/cmath.sqrt(a_sink)
        old_source = source
        old_sink = sink
        source = a_sink * inv_sigma2_sink * mat[ts].H * old_sink
        sink = a_source * inv_sigma2_source * mat[ts] * old_source
        # source = a_source * inv_sigma2_source * mat[ts] * old_source
        source /= cmath.sqrt((source.H * source)[0,0])
        sink /= cmath.sqrt((sink.H * sink)[0,0])
        # print('Iter '+str(n)+' '+ str((old_source.H * source)[0,0])+' '+str((old_sink.H * sink)[0,0])+' '+str((source.H * sink)[0,0]))
    return source, sink

# def find_opt_sn(source_g, sink_g, corr, mat, ts):
#     mat_size, temp = mat[0].shape
#     n_configs, n_types, t_size = corr.shape
#     n_iter = 20
#     source = source_g
#     sink = sink_g
#     c_av = mat[ts]
#     for x in range(mat_size):
#         for y in range(mat_size):
#             if(x == y):
#                 c_av[x,y] = c_av[x,y].real
#             elif(x < y):
#                 re = (c_av[x,y].real+c_av[y,x].real)/2.0
#                 im = (c_av[x,y].imag-c_av[y,x].imag)/2.0
#                 c_av[x,y] = re+im*1.0j
#                 c_av[y,x] = re-im*1.0j
#     for t in range(n_iter):
#         sigma2_source = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
#         sigma2_sink = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
#         for x in range(n_configs):
#             c = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
#             for y in range(n_types):
#                 c[math.floor(y/mat_size), y%mat_size] = corr[x, y, ts]
#             for n in range(mat_size):
#                 for m in range(mat_size):
#                     if(n == m):
#                         c[n,m] = c[n,m].real
#                     elif(n < m):
#                         re = (c[n,m].real+c[m,n].real)/2.0
#                         im = (c[n,m].imag-c[m,n].imag)/2.0
#                         c[n,m] = re+im*1.0j
#                         c[m,n] = re-im*1.0j
#             sigma2_source += c * source * source.H * c.H
#             sigma2_sink += c.H * sink * sink.H * c
#         sigma2_source /= n_configs
#         sigma2_sink /= n_configs
#         inv_sigma2_source = np.linalg.inv(sigma2_source)
#         inv_sigma2_sink = np.linalg.inv(sigma2_sink)
#         a_source = (source.H * c_av.H * inv_sigma2_source * inv_sigma2_source * c_av * source)[0,0]
#         a_sink = (sink.H * c_av * inv_sigma2_sink * inv_sigma2_sink * c_av.H * sink)[0,0]
#         a_source = 1.0/cmath.sqrt(a_source)
#         a_sink = 1.0/cmath.sqrt(a_sink)
#         old_source = source
#         old_sink = sink
#         source = a_sink * inv_sigma2_sink * c_av.H * old_sink
#         sink = a_source * inv_sigma2_source * c_av * old_source
#         # source = a_source * inv_sigma2_source * c_av * old_source
#         source /= cmath.sqrt((source.H * source)[0,0])
#         sink /= cmath.sqrt((sink.H * sink)[0,0])
#         print('Iter '+str(t)+' '+ str((old_source.H * source)[0,0])+' '+str((old_sink.H * sink)[0,0])+' '+str((source.H * sink)[0,0]))
    # return source, sink


def compute_corr(source, sink, mat):
    t_size = len(mat)
    corr = np.zeros(t_size, dtype=np.complex128)
    for x in range(t_size):
        corr[x] = (sink.H * mat[x] * source)[0,0]
    return corr


def find_gnd_masses(corr):
    n_configs, n_types, t_size = corr.shape
    n_boot = 100
    masses = np.zeros((n_boot, 4, t_size))
    # av_corr = average(corr)
    # temp_mat = make_matrices(av_corr)
    # guess = np.matrix([[0],[1],[0],[0],[0]])
    # source_sn, sink_sn = find_opt_sn(guess, guess, corr, temp_mat, 20)
    for x in range(0, n_boot):
        if(x < n_boot-1):
            print('Performing bootstrap... ' + str(math.floor(x/n_boot*100)) + '%', end='\r')
        else:
            print('Performing bootstrap... \033[1;32msuccess\033[0m')
        temp_corr = np.zeros((n_configs, n_types, t_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_corr[y] = corr[r]
        av_corr = average(temp_corr)
        masses[x, 0] = eff_masses(av_corr[6])
        temp_mat = make_matrices(av_corr)
        temp_evals, temp_evecs = find_eigsys(temp_mat)
        masses[x, 1] = eff_masses(temp_evals[0])
        ts = 20
        source = temp_evecs[0][ts]
        source *= 1.0/cmath.sqrt((source.H * source)[0,0])
        sink = find_sink(temp_corr, source, temp_mat, ts)
        temp_corr2 = compute_corr(source, sink, temp_mat)
        masses[x, 2] = eff_masses(temp_corr2)
        guess = np.matrix([[0],[1],[0],[0],[0]])
        source_sn, sink_sn = find_opt_sn(guess, guess, temp_corr, temp_mat, ts)
        # print(str((source_sn.H * guess)[0,0])+' '+str((sink_sn.H * guess)[0,0]))
        temp_corr3 = compute_corr(source_sn, source_sn, temp_mat)
        masses[x, 3] = eff_masses(temp_corr3)
    av_masses = np.zeros((4, t_size))
    for x in range(0, n_boot):
        av_masses += masses[x]
    av_masses /= float(n_boot)
    err_masses = np.zeros((4, t_size))
    for x in range(0, n_boot):
        err_masses += (masses[x]-av_masses)**2
    err_masses = (1./float(n_boot)*err_masses)**(1/2)
    return av_masses, err_masses




pion020 = read_file(FILE_NAMES)
# proton = read_file('test.SP')

mass, err = find_gnd_masses(pion020)
xlab = np.arange(0, 128, 1)
xerl = np.zeros(128)
# plt.style.use('classic')
# l0 = plt.plot(xlab, mass[0], color='r', marker='o')
# rb0 = plt.errorbar(xlab, mass[0], xerr=xerl, yerr=err[0], ecolor='r', fmt=None)
lines1 = plt.errorbar(xlab, mass[1], xerr=xerl, yerr=err[1], color='b')
lines2 = plt.errorbar(xlab, mass[2], xerr=xerl, yerr=err[2], color='g')
# l0 = plt.plot(xlab, mass[3], color='k', marker='s')
# eb3 = plt.errorbar(xlab, mass[3], xerr=xerl, yerr=err[3], ecolor='k', fmt=None)
# lines = plt.errorbar(xlab, mass[0], xerr=xerl, yerr=err[0], fmt='o', color='r', ecolor='g')
# lines = plt.plot(xlab, mass[0], "bs")
# plt.setp(lines, color='r', linewidth=2.0)
plt.show()