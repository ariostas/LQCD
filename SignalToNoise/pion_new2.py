import numpy as np
import scipy.optimize as opt
import math
import cmath
import random
import matplotlib.pyplot as plt
# import seaborn as sns

random.seed()
FILE_NAMES = ["" for x in range(25)]
for x in range(1, 6):
    for y in range(1, 6):
        FILE_NAMES[(x-1)*5+(y-1)] = '/home/arios/Documents/LQCDConfigs/concatenated_new/proton_SrcDG' + str(x) + '_SnkDG' + str(y) + '_Interp4.dat'


def read_file(names):
    f_counter = 0
    n_files = len(names)
    data = np.zeros((1, n_files, 1), dtype=np.complex128)
    for file in names:
        if(f_counter < n_files-1):
            print('Reading correlator data... ' + str(math.floor(f_counter/n_files*100)) + '%', end='\r')
        else:
            print('Reading correlator data... \033[1;32mdone\033[0m')
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
    n_configs = len(corr)
    av = corr[0]
    for x in range(1, n_configs):
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


def construct_matrices(corr):
    n_configs, n_types, t_size = corr.shape
    mat_size = math.sqrt(n_types)
    if(mat_size%1 != 0):
        print("Matrix size error")
    matrices = [[np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128)) for n in range(t_size)] for m in range(n_configs)]
    matrices = np.array(matrices)
    for x in range(n_configs):
        if(x < n_configs-1):
            print('Constructing matrices... ' + str(math.floor(x/n_configs*100)) + '%', end='\r')
        else:
            print('Constructing matrices... \033[1;32mdone\033[0m')
        for y in range(t_size):
            temp_arr = np.zeros((mat_size, mat_size), dtype=np.complex128)
            for z in range(n_types):
                temp_arr[math.floor(z/mat_size), z%mat_size] = corr[x, z, y]
            matrices[x,y] = np.matrix(temp_arr)
    return matrices

def make_hermitian(mat):
    n_configs, t_size, mat_size, temp = mat.shape
    temp_mat = mat
    for x in range(n_configs):
        if(x < n_configs-1):
            print('Making matrices hermitian... ' + str(math.floor(x/n_configs*100)) + '%', end='\r')
        else:
            print('Making matrices hermitian... \033[1;32mdone\033[0m')
        for y in range(t_size):
            for n in range(mat_size):
                for m in range(mat_size):
                    if(n == m):
                        temp_mat[x,y][n,m] = mat[x,y][n,m].real
                    elif(n > m):
                        re = (mat[x,y][n,m].real + mat[x,y][m,n].real)/2.0
                        im = (mat[x,y][n,m].imag - mat[x,y][m,n].imag)/2.0
                        temp_mat[x,y][n,m] = re + im*1.0j
                        temp_mat[x,y][m,n] = re - im*1.0j
    return temp_mat


def find_eigsys(mat):
    t_size, mat_size, temp = mat.shape
    evals = np.zeros((mat_size, t_size), dtype=np.complex128)
    evecs = [[np.matrix(np.zeros((mat_size, 1)), dtype=np.complex128) for n in range(t_size)] for m in range(mat_size)]
    evecs = np.array(evecs)
    init = np.matrix(mat[0])
    init2 = np.matrix(mat[3])
    q1 = np.linalg.cholesky(init)
    q2 = np.linalg.cholesky(init2)
    invq1 = np.linalg.inv(q1)
    invq2 = np.linalg.inv(q2)
    for x in range(t_size):
        temp_evals, temp_evecs = np.linalg.eig(invq1 * np.matrix(mat[x]) * invq1.H)
        temp_evals2, temp_evecs2 = np.linalg.eig(invq2 * np.matrix(mat[x]) * invq2.H)
        temp_evecs = invq1.H * temp_evecs
        temp_evecs2 = invq2.H * temp_evecs2
        order = np.argsort(temp_evals)
        order2 = np.argsort(temp_evals2)
        for y in range(mat_size):
            evals[y, x] = temp_evals[order[mat_size-1-y]]
            evecs[y, x] = temp_evecs2[:,order2[mat_size-1-y]]
    return evals, evecs


def find_sink(source, mat, ts):
    n_configs, t_size, mat_size, temp = mat.shape
    av_mat = average(mat)
    sink = np.matrix(np.zeros((mat_size, 1), dtype=np.complex128))
    sigma2 = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
    for x in range(n_configs):
        sigma2 += np.matrix(mat[x, ts]) * np.matrix(source) * np.matrix(source).H * np.matrix(mat[x, ts]).H
    sigma2 /= n_configs
    inv_sigma2 = np.linalg.inv(sigma2)
    a = (np.matrix(source).H * np.matrix(av_mat[ts]).H * inv_sigma2 * inv_sigma2 * np.matrix(av_mat[ts]) * np.matrix(source))[0,0]
    a = 1.0/cmath.sqrt(a)
    sink = a * inv_sigma2 * np.matrix(av_mat[ts]) * np.matrix(source)
    sink /= cmath.sqrt((sink.H * sink)[0,0])
    return sink


def find_opt_sn(source_g, sink_g, mat, ts):
    n_configs, t_size, mat_size, temp = mat.shape
    av_mat = average(mat)
    n_iter = 20
    source = np.matrix(source_g)
    sink = np.matrix(sink_g)
    for n in range(n_iter):
        sigma2_source = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
        sigma2_sink = np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128))
        for x in range(n_configs):
            sigma2_source += np.matrix(mat[x, ts]) * source * source.H * np.matrix(mat[x, ts]).H
            sigma2_sink += np.matrix(mat[x, ts]).H * sink * sink.H * np.matrix(mat[x, ts])
        sigma2_source /= float(n_configs)
        sigma2_sink /= float(n_configs)
        inv_sigma2_source = np.linalg.inv(sigma2_source)
        inv_sigma2_sink = np.linalg.inv(sigma2_sink)
        a_source = (source.H * np.matrix(av_mat[ts]).H * inv_sigma2_source * inv_sigma2_source * np.matrix(av_mat[ts]) * source)[0,0]
        a_sink = (sink.H * np.matrix(av_mat[ts]) * inv_sigma2_sink * inv_sigma2_sink * np.matrix(av_mat[ts]).H * sink)[0,0]
        a_source = 1.0/cmath.sqrt(a_source)
        a_sink = 1.0/cmath.sqrt(a_sink)
        old_source = source
        old_sink = sink
        source = a_sink * inv_sigma2_sink * np.matrix(av_mat[ts]).H * old_sink
        sink = a_source * inv_sigma2_source * np.matrix(av_mat[ts]) * old_source
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
        corr[x] = (np.matrix(sink).H * np.matrix(mat[x]) * np.matrix(source))[0,0]
    return corr


def find_gnd_masses(corr, mat):
    n_configs, n_types, t_size = corr.shape
    n_configs, t_size, mat_size, temp = mat.shape
    n_boot = 50
    masses = np.zeros((n_boot, 4, t_size))
    # av_corr = average(corr)
    # temp_mat = make_matrices(av_corr)
    # guess = np.matrix([[0],[1],[0],[0],[0]])
    # source_sn, sink_sn = find_opt_sn(guess, guess, corr, temp_mat, 20)
    for x in range(0, n_boot):
        if(x < n_boot-1):
            print('Performing bootstrap... ' + str(math.floor(x/n_boot*100)) + '%', end='\r')
        else:
            print('Performing bootstrap... \033[1;32mdone\033[0m')
        temp_corr = np.zeros((n_configs, n_types, t_size), dtype=np.complex128)
        temp_mat = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_corr[y] = corr[r]
            temp_mat[y] = mat[r]
        av_corr = average(temp_corr)
        masses[x, 0] = eff_masses(av_corr[6])
        av_mat = average(temp_mat)
        temp_evals, temp_evecs = find_eigsys(av_mat)
        masses[x, 1] = eff_masses(temp_evals[0])
        ts = 20
        source = temp_evecs[0, ts]
        source *= 1.0/cmath.sqrt((np.matrix(source).H * np.matrix(source))[0,0])
        sink = find_sink(source, temp_mat, ts)
        temp_corr2 = compute_corr(source, sink, av_mat)
        masses[x, 2] = eff_masses(temp_corr2)
        guess = np.matrix([[0],[1],[0],[0],[0]])
        source_sn, sink_sn = find_opt_sn(guess, guess, temp_mat, ts)
        # print(str((source_sn.H * guess)[0,0])+' '+str((sink_sn.H * guess)[0,0]))
        temp_corr3 = compute_corr(source_sn, sink_sn, av_mat)
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

# sns.set()
# sns.set_context('paper')
# sns.set_style('ticks')
plt.figure(figsize=(8, 2))


pion020_corr = read_file(FILE_NAMES)
pion020_mat = construct_matrices(pion020_corr)
pion020_mat = make_hermitian(pion020_mat)

mass, err = find_gnd_masses(pion020_corr, pion020_mat)
xlab = np.arange(0, 64, 1)
xerl = np.zeros(64)
# plt.style.use('classic')
# l0 = plt.plot(xlab, np.resize(mass[0], 64), color='r', marker='o')
# eb0 = plt.errorbar(xlab, np.resize(mass[0], 64), xerr=xerl, yerr=np.resize(err[0], 64), ecolor='r', fmt=None)
# l1 = plt.plot(xlab, np.resize(mass[1], 64), color='b', marker='o', linewidth=0, label='Variational')
eb1 = plt.errorbar(xlab, np.resize(mass[1], 64), yerr=np.resize(err[1], 64), color='b', ecolor='b', fmt='^', capsize=2, label='Variational')
# l2 = plt.plot(xlab, np.resize(mass[2], 64), color='g', marker='s', linewidth=0, label='Variational source + s/n sink')
eb2 = plt.errorbar(xlab, np.resize(mass[2], 64), yerr=np.resize(err[2], 64), color='g', ecolor='g', fmt='s', capsize=2, label='Variational source + s/n sink')
# l0 = plt.plot(xlab, np.resize(mass[3], 64), color='k', marker='t')
eb3 = plt.errorbar(xlab, np.resize(mass[3], 64), yerr=np.resize(err[3], 64), color='r', ecolor='r', fmt='o', capsize=2, label='s/n source and sink')
# lines = plt.errorbar(xlab, mass[0], xerr=xerl, yerr=err[0], fmt='o', color='r', ecolor='g')
# lines = plt.plot(xlab, mass[0], "bs")
# plt.setp(l1, linewidth=0)
# plt.setp(l2, color='g', linewidth=0)
plt.ylabel('$m_{eff}$', fontsize=16)
plt.xlabel('$ \Delta t $', fontsize=16)
plt.axis([-1, 65, 0.1, 0.35])
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode='expand', borderaxespad=0.)
plt.show()
