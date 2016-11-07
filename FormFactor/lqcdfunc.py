#!/usr/bin/python3
import numpy as np
import math
import cmath
import random
import matplotlib.pyplot as plt
import fits
import multiprocessing

threeptfn_file = 'bar3ptfn/{0}_cur3ptfn_{1}_i{2}_g{3}_qx{4}_qy{5}_qz{6}_pfx{7}_pfy{8}_pfz{9}.{10}.{11}.sh_{12}_sh_{13}.SS'
twoptfn_0_file = 'hadspec/{0}.D{1}.{2}.{3}.sh_{4}_sh_{5}.SS'
twoptfn_other_file = 'hadspec/{0}_px{1}_py{2}_pz{3}.D{4}.{5}.{6}.sh_{7}_sh_{8}.SS'

random.seed()


def read_file(data_type='2ptfn', m=-8999, sources=['DG1_1'], sinks=['DG1_1'], pf=(0,0,0), q=(1,0,0), g=0, i=0, seqsource='NUCL_D_UNPOL', t_sink=16, current='nonlocal', had='proton', file_prefix = '/home/arios/Documents/LQCDConfigs/wil_16_60_aniso/'):
    n_sources = len(sources)
    n_sinks = len(sinks)
    data = np.zeros((1, n_sources*n_sinks, 1), dtype=np.complex128)
    for x in range(n_sources):
        for y in range(n_sinks):
            if(x*n_sinks + y < n_sources*n_sinks-1):
                print('Reading data... ' + str(math.floor((x*n_sinks + y)/(n_sources*n_sinks-1)*100)) + '%', end='\r')
            else:
                print('Reading data... \033[1;32mdone\033[0m')
            filename = ''
            if(data_type == '2ptfn' and pf == (0,0,0)):
                filename = file_prefix+twoptfn_0_file.format(had, m, sources[x], sinks[y], x, y)
            elif(data_type == '2ptfn'):
                filename = file_prefix+twoptfn_other_file.format(had, pf[0], pf[1], pf[2], m, sources[x], sinks[y], x, y)
            elif(data_type == '3ptfn'):
                filename = file_prefix+threeptfn_file.format(current, seqsource, i, g, q[0], q[1], q[2], pf[0], pf[1], pf[2], sources[x], sinks[y], x, y)
            else:
                print('Error: unknown data type')
            f = open(filename, 'r')
            counter = -1
            n_configs = 0
            t_size = 0
            for line in f:
                line = line.strip()
                col = line.split()
                if(counter == -1):
                    n_configs = int(col[0])
                    t_size = int(col[1])
                    if(x == 0 and y == 0):
                        data = np.zeros((n_configs, n_sources*n_sinks, t_size), dtype=np.complex128)
                else:
                    data[math.floor(counter/t_size), x*n_sinks + y, counter%t_size] = float(col[1]) + float(col[2])*1.0j
                counter += 1
            if(counter != n_configs*t_size):
                print('Error: File seems to be corrupted')
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
            print('Constructing matrices... ' + str(math.floor(x/(n_configs-1)*100)) + '%', end='\r')
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
            print('Making matrices hermitian... ' + str(math.floor(x/(n_configs-1)*100)) + '%', end='\r')
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


def find_eigsys(mat, t0):
    t_size, mat_size, temp = mat.shape
    evals = np.zeros((mat_size, t_size), dtype=np.complex128)
    evecs = [[np.matrix(np.zeros((mat_size, 1)), dtype=np.complex128) for n in range(t_size)] for m in range(mat_size)]
    evecs = np.array(evecs)
    init = np.matrix(mat[t0])
    init2 = np.matrix(mat[t0])
    init = np.linalg.inv(init)
    init2 = np.linalg.inv(init2)
    for x in range(t_size):
        temp_evals, temp_evecs = np.linalg.eig(init * np.matrix(mat[x]))
        temp_evals2, temp_evecs2 = np.linalg.eig(init2 * np.matrix(mat[x]))
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


def find_opt_sn(mat, ts):
    n_configs, t_size, mat_size, temp = mat.shape
    av_mat = average(mat)
    n_iter = 20
    n_tries = 2
    converged = False
    for tries in range(n_tries):
        s1 = np.random.random()+0.001
        s2 = np.random.random()+0.001
        s3 = np.random.random()+0.001
        source = np.matrix([[s1], [s2], [s3], [s1-s3], [s2+s3]], dtype=np.complex128)
        # source = np.matrix([[1], [0], [0], [0], [0]], dtype=np.complex128)
        source *= 1.0/cmath.sqrt((source.H * source)[0,0])
        sink = source
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
            # print('Try '+str(tries)+' Iter '+str(n)+' '+ str((old_source.H * source)[0,0])+' '+str((old_sink.H * sink)[0,0])+' '+str((source.H * sink)[0,0]))
            # print('Iter'+str(n)+'  '+str(source)+'\n'+str(sink))
        if(1 - (old_source.H * source)[0,0].real < 0.01):
            converged = True
            break
    if(not converged):
        print('S/N process did not converge                     ')
        return np.matrix([[0],[0],[0]]), np.matrix([[0],[0],[0]])
    return source, sink


def compute_corr(source, sink, mat):
    t_size = len(mat)
    corr = np.zeros(t_size, dtype=np.complex128)
    for x in range(t_size):
        corr[x] = (np.matrix(sink).H * np.matrix(mat[x]) * np.matrix(source))[0,0]
    return corr


def compute_fullcorr(source, sink, mat):
    n_configs, t_size, mat_size, temp = mat.shape
    corr = np.zeros((n_configs, t_size), dtype=np.complex128)
    for x in range(n_configs):
        corr[x] = compute_corr(source, sink, mat[x])
    return corr


def find_gnd_masses(twoptfn, t0, t1, n_boot = 200):
    n_configs, t_size, mat_size, temp = twoptfn.shape
    n_failed = 0
    masses = np.zeros((n_boot, 8, t_size))
    for x in range(0, n_boot):
        if(x < n_boot-1):
            print('Performing mass bootstrapping... ' + str(math.floor(x/(n_boot-1)*100)) + '%', end='\r')
        else:
            print('Performing mass bootstrapping... \033[1;32mdone\033[0m')
        temp_2ptfn = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_2ptfn[y] = twoptfn[r]
        av_2ptfn = average(temp_2ptfn)
        ##############################
        #### Smeared src and snk #####
        ##############################
        smeared_src1 = np.matrix([[1],[0],[0],[0],[0]])
        smeared_src2 = np.matrix([[0],[1],[0],[0],[0]])
        smeared_src3 = np.matrix([[0],[0],[1],[0],[0]])
        smeared_src4 = np.matrix([[0],[0],[0],[1],[0]])
        smeared_src5 = np.matrix([[0],[0],[0],[0],[1]])
        smeared_2ptfn_corr1 = compute_corr(smeared_src1, smeared_src1, av_2ptfn)
        masses[x, 0] = eff_masses(smeared_2ptfn_corr1)
        smeared_2ptfn_corr2 = compute_corr(smeared_src2, smeared_src2, av_2ptfn)
        masses[x, 1] = eff_masses(smeared_2ptfn_corr2)
        smeared_2ptfn_corr3 = compute_corr(smeared_src3, smeared_src3, av_2ptfn)
        masses[x, 2] = eff_masses(smeared_2ptfn_corr3)
        smeared_2ptfn_corr4 = compute_corr(smeared_src4, smeared_src4, av_2ptfn)
        masses[x, 3] = eff_masses(smeared_2ptfn_corr4)
        smeared_2ptfn_corr5 = compute_corr(smeared_src5, smeared_src5, av_2ptfn)
        masses[x, 4] = eff_masses(smeared_2ptfn_corr5)
        ##############################
        ## Variational src and snk ###
        ##############################
        temp_evals, temp_evecs = find_eigsys(av_2ptfn, t0)
        var_src = temp_evecs[0, t1]
        var_src *= 1.0/cmath.sqrt((np.matrix(var_src).H * np.matrix(var_src))[0,0])
        var_2ptfn_corr = compute_corr(var_src, var_src, av_2ptfn)
        masses[x, 5] = eff_masses(var_2ptfn_corr)
        # print(temp_evals[:,t1])
        # print(var_src)
        # print('var: ',var_src)
        ##############################
        #### Var src and S/N snk #####
        ##############################
        varsn_src = var_src
        varsn_snk = find_sink(varsn_src, temp_2ptfn, t1)
        varsn_2ptfn_corr = compute_corr(varsn_src, varsn_snk, av_2ptfn)
        masses[x, 6] = eff_masses(varsn_2ptfn_corr)
        ##############################
        ###### S/N src and snk #######
        ##############################
        sn_src, sn_snk = find_opt_sn(temp_2ptfn, t1)
        if(sn_src[0] == 0 and sn_src[1] == 0 and sn_src[2] == 0):
            n_failed += 1
        else:
            sn_2ptfn_corr = compute_corr(sn_src, sn_snk, av_2ptfn)
            masses[x, 7] = eff_masses(sn_2ptfn_corr)
            # print(sn_src)
            # print('sn: ',sn_src)
        ##############################
    av_masses = np.zeros((8, t_size))
    for x in range(0, n_boot):
        av_masses += masses[x]
    av_masses /= float(n_boot)
    av_masses[7] *= float(n_boot)/float(n_boot-n_failed)
    err_masses = np.zeros((8, t_size))
    for x in range(0, n_boot):
        for y in range(8):
            if(masses[x][y][0] != 0 or y != 7):
                err_masses[y] += (masses[x][y]-av_masses[y])**2
    err_masses = (1./float(n_boot)*err_masses)**(1/2)
    err_masses[7] *= (float(n_boot)/float(n_boot-n_failed))**(1/2)
    return av_masses, err_masses


def find_gnd_massesandcorrs(twoptfn, t0, t1):
    n_configs, t_size, mat_size, temp = twoptfn.shape
    n_boot = 200
    n_failed = 0
    masses = np.zeros((n_boot, 8, t_size))
    corrs = np.zeros((n_boot, 8, n_configs, t_size))
    for x in range(0, n_boot):
        if(x < n_boot-1):
            print('Performing mass bootstrapping... ' + str(math.floor(x/(n_boot-1)*100)) + '%', end='\r')
        else:
            print('Performing mass bootstrapping... \033[1;32mdone\033[0m')
        temp_2ptfn = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_2ptfn[y] = twoptfn[r]
        av_2ptfn = average(temp_2ptfn)
        ##############################
        #### Smeared src and snk #####
        ##############################
        smeared_src1 = np.matrix([[1],[0],[0],[0],[0]])
        smeared_src2 = np.matrix([[0],[1],[0],[0],[0]])
        smeared_src3 = np.matrix([[0],[0],[1],[0],[0]])
        smeared_src4 = np.matrix([[0],[0],[0],[1],[0]])
        smeared_src5 = np.matrix([[0],[0],[0],[0],[1]])
        corrs[x, 0] = compute_fullcorr(smeared_src1, smeared_src1, temp_2ptfn)
        corrs[x, 1] = compute_fullcorr(smeared_src2, smeared_src2, temp_2ptfn)
        corrs[x, 2] = compute_fullcorr(smeared_src3, smeared_src3, temp_2ptfn)
        corrs[x, 3] = compute_fullcorr(smeared_src4, smeared_src4, temp_2ptfn)
        corrs[x, 4] = compute_fullcorr(smeared_src5, smeared_src5, temp_2ptfn)
        smeared_2ptfn_corr1 = compute_corr(smeared_src1, smeared_src1, av_2ptfn)
        masses[x, 0] = eff_masses(smeared_2ptfn_corr1)
        smeared_2ptfn_corr2 = compute_corr(smeared_src2, smeared_src2, av_2ptfn)
        masses[x, 1] = eff_masses(smeared_2ptfn_corr2)
        smeared_2ptfn_corr3 = compute_corr(smeared_src3, smeared_src3, av_2ptfn)
        masses[x, 2] = eff_masses(smeared_2ptfn_corr3)
        smeared_2ptfn_corr4 = compute_corr(smeared_src4, smeared_src4, av_2ptfn)
        masses[x, 3] = eff_masses(smeared_2ptfn_corr4)
        smeared_2ptfn_corr5 = compute_corr(smeared_src5, smeared_src5, av_2ptfn)
        masses[x, 4] = eff_masses(smeared_2ptfn_corr5)
        ##############################
        ## Variational src and snk ###
        ##############################
        temp_evals, temp_evecs = find_eigsys(av_2ptfn, t0)
        var_src = temp_evecs[0, t1]
        var_src *= 1.0/cmath.sqrt((np.matrix(var_src).H * np.matrix(var_src))[0,0])
        corrs[x, 5] = compute_fullcorr(var_src, var_src, temp_2ptfn)
        var_2ptfn_corr = compute_corr(var_src, var_src, av_2ptfn)
        masses[x, 5] = eff_masses(var_2ptfn_corr)
        # print(temp_evals[:,t1])
        print('var: ',var_src)
        ##############################
        #### Var src and S/N snk #####
        ##############################
        varsn_src = var_src
        varsn_snk = find_sink(varsn_src, temp_2ptfn, t1)
        corrs[x, 6] = compute_fullcorr(varsn_src, varsn_snk, temp_2ptfn)
        varsn_2ptfn_corr = compute_corr(varsn_src, varsn_snk, av_2ptfn)
        masses[x, 6] = eff_masses(varsn_2ptfn_corr)
        ##############################
        ###### S/N src and snk #######
        ##############################
        guess = np.matrix([[0],[0],[1]])
        sn_src, sn_snk = find_opt_sn(guess, guess, temp_2ptfn, t1)
        if(sn_src[0] == 0 and sn_src[1] == 0 and sn_src[2] == 0):
            n_failed += 1
        else:
            corrs[x, 7] = compute_fullcorr(sn_src, sn_snk, temp_2ptfn)
            sn_2ptfn_corr = compute_corr(sn_src, sn_snk, av_2ptfn)
            masses[x, 7] = eff_masses(sn_2ptfn_corr)
            print('sn: ',sn_src)
        ##############################
    av_masses = np.zeros((8, t_size))
    av_corrs = np.zeros((8, n_configs, t_size))
    for x in range(0, n_boot):
        av_masses += masses[x]
        av_corrs += corrs[x]
    av_masses /= float(n_boot)
    av_masses[7] *= float(n_boot)/float(n_boot-n_failed)
    av_corrs /= float(n_boot)
    av_corrs[7] *= float(n_boot)/float(n_boot-n_failed)
    err_masses = np.zeros((8, t_size))
    err_corrs = np.zeros((8, n_configs, t_size))
    for x in range(0, n_boot):
        for y in range(8):
            if(corrs[x][y][0] != 0 or y != 7):
                err_masses[y] += (masses[x][y]-av_masses[y])**2
                err_corrs[y] += (corrs[x][y]-av_corrs[y])**2
    err_masses = (1./float(n_boot)*err_masses)**(1/2)
    err_masses[7] *= (float(n_boot)/float(n_boot-n_failed))**(1/2)
    err_corrs = (1./float(n_boot)*err_corrs)**(1/2)
    err_corrs[7] *= (float(n_boot)/float(n_boot-n_failed))**(1/2)
    return av_masses, err_masses, av_corrs, err_corrs


def compute_ff_real(threeptfn, twoptfn_src, twoptfn_snk, t_sink):
    ff = np.zeros(t_sink+1)
    for tau in range(t_sink+1):
        r = twoptfn_snk[t_sink].real/twoptfn_src[t_sink].real
        r *= twoptfn_snk[tau].real/twoptfn_src[tau].real
        r *= twoptfn_src[t_sink-tau].real/twoptfn_snk[t_sink-tau].real
        r = math.sqrt(abs(r))
        r *= threeptfn[tau].real/twoptfn_snk[t_sink].real
        ff[tau] = r
    return ff


def compute_ff_complex(threeptfn, twoptfn_src, twoptfn_snk, t_sink):
    ff = np.zeros(t_sink+1)
    for tau in range(t_sink+1):
        r = twoptfn_snk[t_sink]/twoptfn_src[t_sink]
        r *= twoptfn_snk[tau]/twoptfn_src[tau]
        r *= twoptfn_src[t_sink-tau]/twoptfn_snk[t_sink-tau]
        r = cmath.sqrt(r)
        r *= threeptfn[tau]/twoptfn_snk[t_sink]
        ff[tau] = r.real
    return ff


def find_ff(threeptfn, twoptfn_src, twoptfn_snk, t_sink, t0, t1):
    n_configs, t_size, mat_size, temp = threeptfn.shape
    n_boot = 100
    n_failed = 0
    all_ff = np.zeros((n_boot, 8, t_sink+1))
    for x in range(0, n_boot):
        if(x < n_boot-1):
            print('Performing form factor bootstrapping... ' + str(math.floor(x/(n_boot-1)*100)) + '%', end='\r')
        else:
            print('Performing form factor bootstrapping... \033[1;32mdone\033[0m')
        temp_3ptfn = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        temp_2ptfn_src = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        temp_2ptfn_snk = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_3ptfn[y] = threeptfn[r]
            temp_2ptfn_src[y] = twoptfn_src[r]
            temp_2ptfn_snk[y] = twoptfn_snk[r]
        av_3ptfn = average(temp_3ptfn)
        av_2ptfn_src = average(temp_2ptfn_src)
        av_2ptfn_snk = average(temp_2ptfn_snk)
        ##############################
        #### Smeared src and snk #####
        ##############################
        smeared_src1 = np.matrix([[1],[0],[0],[0],[0]])
        smeared_src2 = np.matrix([[0],[1],[0],[0],[0]])
        smeared_src3 = np.matrix([[0],[0],[1],[0],[0]])
        smeared_src4 = np.matrix([[0],[0],[0],[1],[0]])
        smeared_src5 = np.matrix([[0],[0],[0],[0],[1]])
        smeared_3ptfn_corr1 = compute_corr(smeared_src1, smeared_src1, av_3ptfn)
        smeared_3ptfn_corr2 = compute_corr(smeared_src2, smeared_src2, av_3ptfn)
        smeared_3ptfn_corr3 = compute_corr(smeared_src3, smeared_src3, av_3ptfn)
        smeared_3ptfn_corr4 = compute_corr(smeared_src4, smeared_src4, av_3ptfn)
        smeared_3ptfn_corr5 = compute_corr(smeared_src5, smeared_src5, av_3ptfn)
        smeared_2ptfn_src_corr1 = compute_corr(smeared_src1, smeared_src1, av_2ptfn_src)
        smeared_2ptfn_src_corr2 = compute_corr(smeared_src2, smeared_src2, av_2ptfn_src)
        smeared_2ptfn_src_corr3 = compute_corr(smeared_src3, smeared_src3, av_2ptfn_src)
        smeared_2ptfn_src_corr4 = compute_corr(smeared_src4, smeared_src4, av_2ptfn_src)
        smeared_2ptfn_src_corr5 = compute_corr(smeared_src5, smeared_src5, av_2ptfn_src)
        smeared_2ptfn_snk_corr1 = compute_corr(smeared_src1, smeared_src1, av_2ptfn_snk)
        smeared_2ptfn_snk_corr2 = compute_corr(smeared_src2, smeared_src2, av_2ptfn_snk)
        smeared_2ptfn_snk_corr3 = compute_corr(smeared_src3, smeared_src3, av_2ptfn_snk)
        smeared_2ptfn_snk_corr4 = compute_corr(smeared_src4, smeared_src4, av_2ptfn_snk)
        smeared_2ptfn_snk_corr5 = compute_corr(smeared_src5, smeared_src5, av_2ptfn_snk)
        all_ff[x, 0] = compute_ff_complex(smeared_3ptfn_corr1, smeared_2ptfn_src_corr1, smeared_2ptfn_snk_corr1, t_sink)
        all_ff[x, 1] = compute_ff_complex(smeared_3ptfn_corr2, smeared_2ptfn_src_corr2, smeared_2ptfn_snk_corr2, t_sink)
        all_ff[x, 2] = compute_ff_complex(smeared_3ptfn_corr3, smeared_2ptfn_src_corr3, smeared_2ptfn_snk_corr3, t_sink)
        all_ff[x, 3] = compute_ff_complex(smeared_3ptfn_corr4, smeared_2ptfn_src_corr4, smeared_2ptfn_snk_corr4, t_sink)
        all_ff[x, 4] = compute_ff_complex(smeared_3ptfn_corr5, smeared_2ptfn_src_corr5, smeared_2ptfn_snk_corr5, t_sink)
        ##############################
        ## Variational src and snk ###
        ##############################
        temp_evals, temp_evecs = find_eigsys(av_2ptfn_snk, t0)
        var_src = temp_evecs[0, t1]
        var_src *= 1.0/cmath.sqrt((np.matrix(var_src).H * np.matrix(var_src))[0,0])
        # print(var_src)
        var_3ptfn_corr = compute_corr(var_src, var_src, av_3ptfn)
        var_2ptfn_src_corr = compute_corr(var_src, var_src, av_2ptfn_src)
        var_2ptfn_snk_corr = compute_corr(var_src, var_src, av_2ptfn_snk)
        all_ff[x, 5] = compute_ff_complex(var_3ptfn_corr, var_2ptfn_src_corr, var_2ptfn_snk_corr, t_sink)
        ##############################
        #### Var src and S/N snk #####
        ##############################
        varsn_snk = find_sink(var_src, temp_3ptfn, t1)
        varsn_3ptfn_corr = compute_corr(var_src, varsn_snk, av_3ptfn)
        varsn_2ptfn_src_corr = compute_corr(var_src, varsn_snk, av_2ptfn_src)
        varsn_2ptfn_snk_corr = compute_corr(var_src, varsn_snk, av_2ptfn_snk)
        all_ff[x, 6] = compute_ff_complex(varsn_3ptfn_corr, varsn_2ptfn_src_corr, varsn_2ptfn_snk_corr, t_sink)
        ##############################
        ###### S/N src and snk #######
        ##############################
        sn_src, sn_snk = find_opt_sn(temp_3ptfn, t1)
        if(sn_src[0] == 0 and sn_src[1] == 0 and sn_src[2] == 0):
            n_failed += 1
        else:
            sn_3ptfn_corr = compute_corr(sn_src, sn_snk, av_3ptfn)
            sn_2ptfn_src_corr = compute_corr(sn_src, sn_snk, av_2ptfn_src)
            sn_2ptfn_snk_corr = compute_corr(sn_src, sn_snk, av_2ptfn_snk)
            all_ff[x, 7] = compute_ff_complex(sn_3ptfn_corr, sn_2ptfn_src_corr, sn_2ptfn_snk_corr, t_sink)
        ##############################
    av_ff = np.zeros((8, t_sink+1))
    for x in range(0, n_boot):
        av_ff += all_ff[x]
    av_ff /= float(n_boot)
    av_ff[7] *= float(n_boot)/float(n_boot-n_failed)
    err_ff = np.zeros((8, t_sink+1))
    for x in range(0, n_boot):
        for y in range(8):
            if(all_ff[x][y][0] != 0 or y != 3):
                err_ff[y] += (all_ff[x][y]-av_ff[y])**2
    err_ff = (1./float(n_boot)*err_ff)**(1/2)
    err_ff[7] *= (float(n_boot)/float(n_boot-n_failed))**(1/2)
    return av_ff, err_ff


def find_ff_with2ptfn(threeptfn, twoptfn_src, twoptfn_snk, t_sink, t0, t1):
    n_configs, t_size, mat_size, temp = threeptfn.shape
    n_boot = 100
    n_failed = 0
    all_ff = np.zeros((n_boot, 8, t_sink+1))
    for x in range(0, n_boot):
        if(x < n_boot-1):
            print('Performing form factor bootstrapping... ' + str(math.floor(x/(n_boot-1)*100)) + '%', end='\r')
        else:
            print('Performing form factor bootstrapping... \033[1;32mdone\033[0m')
        temp_3ptfn = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        temp_2ptfn_src = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        temp_2ptfn_snk = np.zeros((n_configs, t_size, mat_size, mat_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_3ptfn[y] = threeptfn[r]
            temp_2ptfn_src[y] = twoptfn_src[r]
            temp_2ptfn_snk[y] = twoptfn_snk[r]
        av_3ptfn = average(temp_3ptfn)
        av_2ptfn_src = average(temp_2ptfn_src)
        av_2ptfn_snk = average(temp_2ptfn_snk)
        ##############################
        #### Smeared src and snk #####
        ##############################
        smeared_src1 = np.matrix([[1],[0],[0],[0],[0]])
        smeared_src2 = np.matrix([[0],[1],[0],[0],[0]])
        smeared_src3 = np.matrix([[0],[0],[1],[0],[0]])
        smeared_src4 = np.matrix([[0],[0],[0],[1],[0]])
        smeared_src5 = np.matrix([[0],[0],[0],[0],[1]])
        smeared_3ptfn_corr1 = compute_corr(smeared_src1, smeared_src1, av_3ptfn)
        smeared_3ptfn_corr2 = compute_corr(smeared_src2, smeared_src2, av_3ptfn)
        smeared_3ptfn_corr3 = compute_corr(smeared_src3, smeared_src3, av_3ptfn)
        smeared_3ptfn_corr4 = compute_corr(smeared_src4, smeared_src4, av_3ptfn)
        smeared_3ptfn_corr5 = compute_corr(smeared_src5, smeared_src5, av_3ptfn)
        smeared_2ptfn_src_corr1 = compute_corr(smeared_src1, smeared_src1, av_2ptfn_src)
        smeared_2ptfn_src_corr2 = compute_corr(smeared_src2, smeared_src2, av_2ptfn_src)
        smeared_2ptfn_src_corr3 = compute_corr(smeared_src3, smeared_src3, av_2ptfn_src)
        smeared_2ptfn_src_corr4 = compute_corr(smeared_src4, smeared_src4, av_2ptfn_src)
        smeared_2ptfn_src_corr5 = compute_corr(smeared_src5, smeared_src5, av_2ptfn_src)
        smeared_2ptfn_snk_corr1 = compute_corr(smeared_src1, smeared_src1, av_2ptfn_snk)
        smeared_2ptfn_snk_corr2 = compute_corr(smeared_src2, smeared_src2, av_2ptfn_snk)
        smeared_2ptfn_snk_corr3 = compute_corr(smeared_src3, smeared_src3, av_2ptfn_snk)
        smeared_2ptfn_snk_corr4 = compute_corr(smeared_src4, smeared_src4, av_2ptfn_snk)
        smeared_2ptfn_snk_corr5 = compute_corr(smeared_src5, smeared_src5, av_2ptfn_snk)
        all_ff[x, 0] = compute_ff_complex(smeared_3ptfn_corr1, smeared_2ptfn_src_corr1, smeared_2ptfn_snk_corr1, t_sink)
        all_ff[x, 1] = compute_ff_complex(smeared_3ptfn_corr2, smeared_2ptfn_src_corr2, smeared_2ptfn_snk_corr2, t_sink)
        all_ff[x, 2] = compute_ff_complex(smeared_3ptfn_corr3, smeared_2ptfn_src_corr3, smeared_2ptfn_snk_corr3, t_sink)
        all_ff[x, 3] = compute_ff_complex(smeared_3ptfn_corr4, smeared_2ptfn_src_corr4, smeared_2ptfn_snk_corr4, t_sink)
        all_ff[x, 4] = compute_ff_complex(smeared_3ptfn_corr5, smeared_2ptfn_src_corr5, smeared_2ptfn_snk_corr5, t_sink)
        ##############################
        ## Variational src and snk ###
        ##############################
        temp_evals, temp_evecs = find_eigsys(av_2ptfn_snk, t0)
        var_src = temp_evecs[0, t1]
        var_src *= 1.0/cmath.sqrt((np.matrix(var_src).H * np.matrix(var_src))[0,0])
        # print(var_src)
        var_3ptfn_corr = compute_corr(var_src, var_src, av_3ptfn)
        var_2ptfn_src_corr = compute_corr(var_src, var_src, av_2ptfn_src)
        var_2ptfn_snk_corr = compute_corr(var_src, var_src, av_2ptfn_snk)
        all_ff[x, 5] = compute_ff_complex(var_3ptfn_corr, var_2ptfn_src_corr, var_2ptfn_snk_corr, t_sink)
        ##############################
        #### Var src and S/N snk #####
        ##############################
        varsn_snk = find_sink(var_src, temp_2ptfn_snk, t1)
        varsn_3ptfn_corr = compute_corr(var_src, varsn_snk, av_3ptfn)
        varsn_2ptfn_src_corr = compute_corr(var_src, varsn_snk, av_2ptfn_src)
        varsn_2ptfn_snk_corr = compute_corr(var_src, varsn_snk, av_2ptfn_snk)
        all_ff[x, 6] = compute_ff_complex(varsn_3ptfn_corr, varsn_2ptfn_src_corr, varsn_2ptfn_snk_corr, t_sink)
        ##############################
        ###### S/N src and snk #######
        ##############################
        sn_src, sn_snk = find_opt_sn(temp_2ptfn_snk, t1)
        if(sn_src[0] == 0 and sn_src[1] == 0 and sn_src[2] == 0):
            n_failed += 1
        else:
            sn_3ptfn_corr = compute_corr(sn_src, sn_snk, av_3ptfn)
            sn_2ptfn_src_corr = compute_corr(sn_src, sn_snk, av_2ptfn_src)
            sn_2ptfn_snk_corr = compute_corr(sn_src, sn_snk, av_2ptfn_snk)
            all_ff[x, 7] = compute_ff_complex(sn_3ptfn_corr, sn_2ptfn_src_corr, sn_2ptfn_snk_corr, t_sink)
        ##############################
    av_ff = np.zeros((8, t_sink+1))
    for x in range(0, n_boot):
        av_ff += all_ff[x]
    av_ff /= float(n_boot)
    av_ff[7] *= float(n_boot)/float(n_boot-n_failed)
    err_ff = np.zeros((8, t_sink+1))
    for x in range(0, n_boot):
        for y in range(8):
            if(all_ff[x][y][0] != 0 or y != 3):
                err_ff[y] += (all_ff[x][y]-av_ff[y])**2
    err_ff = (1./float(n_boot)*err_ff)**(1/2)
    err_ff[7] *= (float(n_boot)/float(n_boot-n_failed))**(1/2)
    return av_ff, err_ff