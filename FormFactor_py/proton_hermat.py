#!/usr/bin/python3
import numpy as np
import math
import cmath
import random
import matplotlib.pyplot as plt
import fits
import multiprocessing

# file_prefix = '/home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450/'
file_prefix = '/home/arios/Documents/LQCDConfigs/wil_16_64_aniso/'
threeptfn_file = 'bar3ptfn/{0}_cur3ptfn_{1}_i{2}_g{3}_qx{4}_qy{5}_qz{6}_pfx{7}_pfy{8}_pfz{9}.{10}.{11}.sh_{12}_sh_{13}.SS'
twoptfn_0_file = 'hadspec/{0}.D{1}.{2}.{3}.sh_{4}_sh_{5}.SS'
twoptfn_other_file = 'hadspec/{0}_px{1}_py{2}_pz{3}.D{4}.{5}.{6}.sh_{7}_sh_{8}.SS'


def read_file(data_type='2ptfn', m=-8999, sources=['DG1_1'], sinks=['DG1_1'], pf=(0,0,0), q=(1,0,0), g=0, i=0, seqsource='NUCL_D_UNPOL', t_sink=16, current='nonlocal', had='proton'):
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


def construct_gen_matrices(corr3, corr2_src, corr2_snk, tau_min, tau_max, t_sink):
    n_configs, n_types, t_size = corr3.shape
    n_blocks = math.sqrt(n_types)
    n_times = tau_max-tau_min+1
    mat_size = math.sqrt(n_types)*(tau_max-tau_min+1)
    if(n_blocks%1 != 0):
        print("Matrix size error")
    matrices = [np.matrix(np.zeros((mat_size, mat_size), dtype=np.complex128)) for n in range(n_configs)]
    matrices = np.array(matrices)
    for x in range(n_configs):
        if(x < n_configs-1):
            print('Constructing matrices... ' + str(math.floor(x/(n_configs-1)*100)) + '%', end='\r')
        else:
            print('Constructing matrices... \033[1;32mdone\033[0m')
        for tau in range(tau_min, tau_max+1):
            for i in range(int(n_blocks)):
                for j in range(int(n_blocks)):
                    temp_entry = corr2_snk[x, i*n_blocks+j, t_sink]/corr2_src[x, i*n_blocks+j, t_sink]
                    temp_entry *= corr2_snk[x, i*n_blocks+j, tau]/corr2_src[x, i*n_blocks+j, tau]
                    temp_entry *= corr2_src[x, i*n_blocks+j, t_sink-tau]/corr2_snk[x, i*n_blocks+j, t_sink-tau]
                    temp_entry = math.sqrt(abs(temp_entry))
                    temp_entry *= abs(corr3[x, i*n_blocks+j, tau])/abs(corr2_snk[x, i*n_blocks+j, t_sink])
                    matrices[x,n_times*i+(tau-tau_min),n_times*j+(tau-tau_min)] = temp_entry
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


def make_hermitian_gen(mat):
    n_configs, mat_size, temp = mat.shape
    temp_mat = mat
    for x in range(n_configs):
        if(x < n_configs-1):
            print('Making matrices hermitian... ' + str(math.floor(x/(n_configs-1)*100)) + '%', end='\r')
        else:
            print('Making matrices hermitian... \033[1;32mdone\033[0m')
        for n in range(mat_size):
            for m in range(mat_size):
                if(n == m):
                    temp_mat[x][n,m] = mat[x][n,m].real
                elif(n > m):
                    re = (mat[x][n,m].real + mat[x][m,n].real)/2.0
                    im = (mat[x][n,m].imag - mat[x][m,n].imag)/2.0
                    temp_mat[x][n,m] = re + im*1.0j
                    temp_mat[x][m,n] = re - im*1.0j
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


def find_gen_src(sink, mat):
    n_configs, mat_ver_size, mat_hor_size = mat.shape
    av_mat = average(mat)
    source = np.matrix(np.zeros((mat_hor_size, 1), dtype=np.complex128))
    sigma2 = np.matrix(np.zeros((mat_hor_size, mat_hor_size), dtype=np.complex128))
    for x in range(n_configs):
        sigma2 += np.matrix(mat[x]).H * np.matrix(sink) * np.matrix(sink).H * np.matrix(mat[x])
    sigma2 /= n_configs
    inv_sigma2 = np.linalg.inv(sigma2)
    a = (np.matrix(sink).H * np.matrix(av_mat) * inv_sigma2 * inv_sigma2 * np.matrix(av_mat).H * np.matrix(sink))[0,0]
    a = 1.0/cmath.sqrt(a)
    source = a * inv_sigma2 * np.matrix(av_mat).H * np.matrix(sink)
    source /= cmath.sqrt((source.H * source)[0,0])
    return source


def find_opt_sn(source_g, sink_g, mat, ts):
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


def find_gen_opt_sn(mat):
    n_configs, mat_ver_size, mat_hor_size = mat.shape
    if(mat_ver_size != mat_hor_size):
        print('Matrix size error')
    av_mat = average(mat)
    evals, evecs = np.linalg.eig(av_mat.real)
    evecs = np.matrix(evecs)
    badevecs_i = []
    n_goodevals = 0
    for i in range(len(evals)):
        ev = evals[i]
        if(ev < 0.5):
            badevecs_i.append(i)
        else:
            n_goodevals += 1
    print('N good evals = ',n_goodevals)
    badevecs = [[evecs[:,i]] for i in badevecs_i]
    n_iter = 20
    n_tries = 3
    converged = False
    for tries in range(n_tries):
        source = np.matrix(np.zeros((mat_hor_size,1), dtype=np.complex128))
        for i in range(mat_hor_size):
            source[i] = np.random.random()+0.001
        source *= 1.0/cmath.sqrt((source.H * source)[0,0])
        sink = np.matrix(np.zeros((mat_ver_size,1), dtype=np.complex128))
        for i in range(mat_ver_size):
            sink[i] = np.random.random()+0.001
        sink *= 1.0/cmath.sqrt((sink.H * sink)[0,0])
        for n in range(n_iter):
            sigma2_source = np.matrix(np.zeros((mat_ver_size, mat_ver_size), dtype=np.complex128))
            sigma2_sink = np.matrix(np.zeros((mat_hor_size, mat_hor_size), dtype=np.complex128))
            for x in range(n_configs):
                sigma2_source += np.matrix(mat[x]) * source * source.H * np.matrix(mat[x]).H
                sigma2_sink += np.matrix(mat[x]).H * sink * sink.H * np.matrix(mat[x])
            sigma2_source /= float(n_configs)
            sigma2_sink /= float(n_configs)
            inv_sigma2_source = np.linalg.inv(sigma2_source)
            inv_sigma2_sink = np.linalg.inv(sigma2_sink)
            a_source = (source.H * np.matrix(av_mat).H * inv_sigma2_source * inv_sigma2_source * np.matrix(av_mat) * source)[0,0]
            a_sink = (sink.H * np.matrix(av_mat) * inv_sigma2_sink * inv_sigma2_sink * np.matrix(av_mat).H * sink)[0,0]
            a_source = 1.0/cmath.sqrt(a_source)
            a_sink = 1.0/cmath.sqrt(a_sink)
            old_source = source
            old_sink = sink
            source = a_sink * inv_sigma2_sink * np.matrix(av_mat).H * old_sink
            sink = a_source * inv_sigma2_source * np.matrix(av_mat) * old_source
            source /= cmath.sqrt((source.H * source)[0,0])
            sink /= cmath.sqrt((sink.H * sink)[0,0])
            for x in range(len(badevecs)):
                evec = np.matrix(badevecs[x][0])
                source = source - (evec.H * source)[0,0]/(evec.H * evec)[0,0] * evec
                sink = sink - (evec.H * sink)[0,0]/(evec.H * evec)[0,0] * evec
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
    return corr.real


def find_gnd_masses(twoptfn, t0, t1):
    n_configs, t_size, mat_size, temp = twoptfn.shape
    n_boot = 200
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
        guess = np.matrix([[0],[0],[1]])
        sn_src, sn_snk = find_opt_sn(guess, guess, temp_2ptfn, t1)
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


def compute_ff(threeptfn, twoptfn_src, twoptfn_snk, t_sink):
    ff = np.zeros(t_sink+1)
    for tau in range(t_sink+1):
        r = twoptfn_snk[t_sink].real/twoptfn_src[t_sink].real
        r *= twoptfn_snk[tau].real/twoptfn_src[tau].real
        r *= twoptfn_src[t_sink-tau].real/twoptfn_snk[t_sink-tau].real
        r = math.sqrt(abs(r))
        r *= threeptfn[tau].real/twoptfn_snk[t_sink].real
        ff[tau] = r
    return ff


def find_ff(threeptfn, twoptfn_src, twoptfn_snk, t_sink, t0, t1):
    n_configs, t_size, mat_size, temp = threeptfn.shape
    n_boot = 50
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
        all_ff[x, 0] = compute_ff(smeared_3ptfn_corr1, smeared_2ptfn_src_corr1, smeared_2ptfn_snk_corr1, t_sink)
        all_ff[x, 1] = compute_ff(smeared_3ptfn_corr2, smeared_2ptfn_src_corr2, smeared_2ptfn_snk_corr2, t_sink)
        all_ff[x, 2] = compute_ff(smeared_3ptfn_corr3, smeared_2ptfn_src_corr3, smeared_2ptfn_snk_corr3, t_sink)
        all_ff[x, 3] = compute_ff(smeared_3ptfn_corr4, smeared_2ptfn_src_corr4, smeared_2ptfn_snk_corr4, t_sink)
        all_ff[x, 4] = compute_ff(smeared_3ptfn_corr5, smeared_2ptfn_src_corr5, smeared_2ptfn_snk_corr5, t_sink)
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
        all_ff[x, 5] = compute_ff(var_3ptfn_corr, var_2ptfn_src_corr, var_2ptfn_snk_corr, t_sink)
        ##############################
        #### Var src and S/N snk #####
        ##############################
        varsn_snk = find_sink(var_src, temp_3ptfn, t1)
        varsn_3ptfn_corr = compute_corr(var_src, varsn_snk, av_3ptfn)
        varsn_2ptfn_src_corr = compute_corr(var_src, varsn_snk, av_2ptfn_src)
        varsn_2ptfn_snk_corr = compute_corr(var_src, varsn_snk, av_2ptfn_snk)
        all_ff[x, 6] = compute_ff(varsn_3ptfn_corr, varsn_2ptfn_src_corr, varsn_2ptfn_snk_corr, t_sink)
        ##############################
        ###### S/N src and snk #######
        ##############################
        guess = var_src
        sn_src, sn_snk = find_opt_sn(guess, guess, temp_3ptfn, t1)
        if(sn_src[0] == 0 and sn_src[1] == 0 and sn_src[2] == 0):
            n_failed += 1
        else:
            sn_3ptfn_corr = compute_corr(sn_src, sn_snk, av_3ptfn)
            sn_2ptfn_src_corr = compute_corr(sn_src, sn_snk, av_2ptfn_src)
            sn_2ptfn_snk_corr = compute_corr(sn_src, sn_snk, av_2ptfn_snk)
            all_ff[x, 7] = compute_ff(sn_3ptfn_corr, sn_2ptfn_src_corr, sn_2ptfn_snk_corr, t_sink)
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


def find_gen_ff(threeptfn, twoptfn_src, twoptfn_snk, ffratio, t_sink, t0, t1):
    n_configs, t_size, mat_size, temp = threeptfn.shape
    n_configs, mat_ver_size, mat_hor_size = ffratio.shape
    n_boot = 50
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
        temp_ffratio = np.zeros((n_configs, mat_ver_size, mat_hor_size), dtype=np.complex128)
        for y in range(0, n_configs):
            r = random.randrange(0, n_configs)
            temp_3ptfn[y] = threeptfn[r]
            temp_2ptfn_src[y] = twoptfn_src[r]
            temp_2ptfn_snk[y] = twoptfn_snk[r]
            temp_ffratio[y] = ffratio[r]
        av_3ptfn = average(temp_3ptfn)
        av_2ptfn_src = average(temp_2ptfn_src)
        av_2ptfn_snk = average(temp_2ptfn_snk)
        av_ffratio = average(temp_ffratio)
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
        all_ff[x, 0] = compute_ff(smeared_3ptfn_corr1, smeared_2ptfn_src_corr1, smeared_2ptfn_snk_corr1, t_sink)
        all_ff[x, 1] = compute_ff(smeared_3ptfn_corr2, smeared_2ptfn_src_corr2, smeared_2ptfn_snk_corr2, t_sink)
        all_ff[x, 2] = compute_ff(smeared_3ptfn_corr3, smeared_2ptfn_src_corr3, smeared_2ptfn_snk_corr3, t_sink)
        all_ff[x, 3] = compute_ff(smeared_3ptfn_corr4, smeared_2ptfn_src_corr4, smeared_2ptfn_snk_corr4, t_sink)
        all_ff[x, 4] = compute_ff(smeared_3ptfn_corr5, smeared_2ptfn_src_corr5, smeared_2ptfn_snk_corr5, t_sink)
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
        all_ff[x, 5] = compute_ff(var_3ptfn_corr, var_2ptfn_src_corr, var_2ptfn_snk_corr, t_sink)
        ##############################
        #### Var snk and S/N src #####
        ##############################
        # varsn_src = find_gen_src(var_src, temp_ffratio)
        # varsn_ff = (np.matrix(var_src).H * np.matrix(av_ffratio) * np.matrix(varsn_src))[0,0].real
        U, s, V = np.linalg.svd(np.matrix(av_ffratio))
        print('svd = ')
        print(s)
        print('svd^2 =')
        print(s**2)
        print('av svd =')
        print(average(s))
        print('av svd^2 =')
        print(average(s**2))
        print('\n')
        # print(varsn_src)
        # print(varsn_ff)
        print("\n\n")
        # all_ff[x, 6] = all_ff[x, 6] + varsn_ff
        all_ff[x, 6] = average(s**2)/float(mat_hor_size)
        ##############################
        ###### S/N src and snk #######
        ##############################
        sn_src, sn_snk = find_gen_opt_sn(temp_ffratio)
        if(sn_src[0] == 0 and sn_src[1] == 0 and sn_src[2] == 0):
            n_failed += 1
        else:
            sn_ff = (np.matrix(sn_snk).H * np.matrix(av_ffratio) * np.matrix(sn_src))[0,0]
            all_ff[x, 7] = all_ff[x, 7] + sn_ff
        ##############################
    av_ff = np.zeros((8, t_sink+1))
    for x in range(0, n_boot):
        av_ff += all_ff[x]
    av_ff /= float(n_boot)
    if(n_failed != n_boot):
        av_ff[7] *= float(n_boot)/float(n_boot-n_failed)
    err_ff = np.zeros((8, t_sink+1))
    for x in range(0, n_boot):
        for y in range(8):
            if(all_ff[x][y][0] != 0 or y != 3):
                err_ff[y] += (all_ff[x][y]-av_ff[y])**2
    err_ff = (1./float(n_boot)*err_ff)**(1/2)
    if(n_failed != n_boot):
        err_ff[7] *= (float(n_boot)/float(n_boot-n_failed))**(1/2)
    return av_ff, err_ff


random.seed()

sources = ['DG0_1', 'DG1_1', 'DG1_1', 'DG2_1', 'DG2_1']
sinks = ['DG0_1', 'DG1_1', 'DG1_1', 'DG2_1', 'DG2_1']

# t_0 = 0 and t_1 = 6 work good
t0 = 2 # time used for generalized eigenvalue problem
t1 = 8 # time at which variational source is picked, and used for S/N optimization
t_sink = 16

# proton_3ptfn_x_corr = read_file(data_type='3ptfn', g=1, sources=sources, sinks=sinks, q=(1,0,0), pf=(0,0,0), current='nonlocal')
# proton_3ptfn_y_corr = read_file(data_type='3ptfn', g=2, sources=sources, sinks=sinks, q=(1,0,0), pf=(0,0,0), current='nonlocal')
# proton_3ptfn_z_corr = read_file(data_type='3ptfn', g=4, sources=sources, sinks=sinks, q=(1,0,0), pf=(0,0,0), current='nonlocal')
proton_3ptfn_t_corr = read_file(data_type='3ptfn', g=8, sources=sources, sinks=sinks, q=(1,0,0), pf=(0,0,0), current='local')
proton_2ptfn_srcp_corr = read_file(data_type='2ptfn', had='proton', pf=(-1,0,0), sources=sources, sinks=sinks)
proton_2ptfn_snkp_corr = read_file(data_type='2ptfn', had='proton', pf=(0,0,0), sources=sources, sinks=sinks)

# proton_3ptfn_x_mat = construct_matrices(proton_3ptfn_x_corr)
# proton_3ptfn_y_mat = construct_matrices(proton_3ptfn_y_corr)
# proton_3ptfn_z_mat = construct_matrices(proton_3ptfn_z_corr)
proton_3ptfn_t_mat = construct_matrices(proton_3ptfn_t_corr)
proton_2ptfn_srcp_mat = construct_matrices(proton_2ptfn_srcp_corr)
proton_2ptfn_snkp_mat = construct_matrices(proton_2ptfn_snkp_corr)
proton_ffratio_mat = construct_gen_matrices(proton_3ptfn_t_corr, proton_2ptfn_srcp_corr, proton_2ptfn_snkp_corr, 8, 10, 16)

proton_ffratio_mat = make_hermitian_gen(proton_ffratio_mat)

# proton_3ptfn_x_mat = make_hermitian(proton_3ptfn_x_mat)
# proton_3ptfn_y_mat = make_hermitian(proton_3ptfn_y_mat)
# proton_3ptfn_z_mat = make_hermitian(proton_3ptfn_z_mat)
# proton_3ptfn_t_mat = make_hermitian(proton_3ptfn_t_mat)
# proton_2ptfn_srcp_mat = make_hermitian(proton_2ptfn_srcp_mat)
# proton_2ptfn_snkp_mat = make_hermitian(proton_2ptfn_snkp_mat)

# proton_ff_x, proton_fferr_x = find_ff(proton_3ptfn_x_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, t_sink, t0, t1)
# proton_ff_y, proton_fferr_y = find_ff(proton_3ptfn_y_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, t_sink, t0, t1)
# proton_ff_z, proton_fferr_z = find_ff(proton_3ptfn_z_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, t_sink, t0, t1)
# proton_ff_t, proton_fferr_t = find_ff(proton_3ptfn_t_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, t_sink, t0, t1)
# rho_mass, proton_masserr = find_gnd_masses(rho_2ptfn_snkp_mat, t0, t1)

proton_ff_t, proton_fferr_t = find_gen_ff(proton_3ptfn_t_mat, proton_2ptfn_srcp_mat, proton_2ptfn_snkp_mat, proton_ffratio_mat, t_sink, t0, t1)

# smeared_src1 = np.matrix([[1],[0],[0],[0],[0]])
# smeared_src2 = np.matrix([[0],[1],[0],[0],[0]])
# smeared_src3 = np.matrix([[0],[0],[1],[0],[0]])
# smeared_src4 = np.matrix([[0],[0],[0],[1],[0]])
# smeared_src5 = np.matrix([[0],[0],[0],[0],[1]])
# corr_src1 = compute_fullcorr(smeared_src1, smeared_src1, rho_2ptfn_snkp_mat)
# corr_src2 = compute_fullcorr(smeared_src2, smeared_src2, rho_2ptfn_snkp_mat)
# corr_src3 = compute_fullcorr(smeared_src3, smeared_src3, rho_2ptfn_snkp_mat)
# corr_src4 = compute_fullcorr(smeared_src4, smeared_src4, rho_2ptfn_snkp_mat)
# corr_src5 = compute_fullcorr(smeared_src5, smeared_src5, rho_2ptfn_snkp_mat)

# av_2ptfn = average(rho_2ptfn_snkp_mat)
# temp_evals, temp_evecs = find_eigsys(av_2ptfn, t0)
# var_src = temp_evecs[0, t1]
# var_src *= 1.0/cmath.sqrt((np.matrix(var_src).H * np.matrix(var_src))[0,0])
# corr_var = compute_fullcorr(var_src, var_src, rho_2ptfn_snkp_mat)

# varsn_src = var_src
# varsn_snk = find_sink(varsn_src, rho_2ptfn_snkp_mat, t1)
# corr_varsn = compute_fullcorr(varsn_src, varsn_snk, rho_2ptfn_snkp_mat)

# guess = np.matrix([[0],[0],[1]])
# sn_src, sn_snk = find_opt_sn(guess, guess, rho_2ptfn_snkp_mat, t1)
# corr_sn = compute_fullcorr(sn_src, sn_snk, rho_2ptfn_snkp_mat)

# masses = np.zeros(8)
# massesstaterr = np.zeros(8)
# massessyserr = np.zeros(8)
# fitrange = np.zeros((8,8))

# results = multiprocessing.Queue()

# threads = []
# threadnames = ['src1', 'src2', 'src3', 'src4', 'src5', 'var', 'varsn', 'sn']
# corrs = [corr_src1, corr_src2, corr_src3, corr_src4, corr_src5, corr_var, corr_varsn, corr_sn]

# for x in range(len(threadnames)):
#     threads.append(multiprocessing.Process(target=fits.scan_fit, args=(corrs[x], [0.0005, 0.5], results, threadnames[x])))

# for x in threads:
#     x.start()

# rho_mass, rho_masserr = find_gnd_masses(rho_2ptfn_snkp_mat, t0, t1)

# for x in threads:
#     x.join()

# for x in range(len(threadnames)):
#     l = results.get()
#     for y in range(len(threadnames)):
#         if(l[0] == threadnames[y]):
#             name, masses[y], massesstaterr[y], massessyserr[y], fitrange[y][0], fitrange[y][1] = l
#             break


# for x in range(8):
#     print(threadnames[x]+': Mass = '+str(masses[x])+' +- '+str(massesstaterr[x])+' (stat) +- '+str(massessyserr[x])+' (sys)')

# masses[5], masseserr[5], fitrange[5][0], fitrange[5][1] = async_result1.get()
# masses[6], masseserr[6], fitrange[6][0], fitrange[6][1] = async_result2.get()
# masses[7], masseserr[7], fitrange[7][0], fitrange[7][1] = async_result3.get()

# masses[0], masseserr[0], fitrange[0][0], fitrange[0][1] = fits.scan_fit(corr_src1, [0.0005, 0.5])
# masses[1], masseserr[1], fitrange[1][0], fitrange[1][1] = fits.scan_fit(corr_src2, [0.0005, 0.5])
# masses[2], masseserr[2], fitrange[2][0], fitrange[2][1] = fits.scan_fit(corr_src3, [0.0005, 0.5])
# masses[3], masseserr[3], fitrange[3][0], fitrange[3][1] = fits.scan_fit(corr_src4, [0.0005, 0.5])
# masses[4], masseserr[4], fitrange[4][0], fitrange[4][1] = fits.scan_fit(corr_src5, [0.0005, 0.5])
# masses[5], masseserr[5], fitrange[5][0], fitrange[5][1] = fits.scan_fit(corr_var, [0.0005, 0.5])
# masses[6], masseserr[6], fitrange[6][0], fitrange[6][1] = fits.scan_fit(corr_varsn, [0.0005, 0.5])
# masses[7], masseserr[7], fitrange[7][0], fitrange[7][1] = fits.scan_fit(corr_sn, [0.0005, 0.5])

# nmin = 7
# nmax= 27
# res, err, chi2, chi2err = fits.fit_correrr(corr, nmin, nmax, [0.0005, 0.5], 'exp')
# print('Mass = '+str(res[1])+' +- '+str(err))
# fits.scan_fit(corr, [0.0005, 0.5])

# proton_ff = [proton_ff_x, proton_ff_y, proton_ff_z, proton_ff_t]
# proton_fferr = [proton_fferr_x, proton_fferr_y, proton_fferr_z, proton_fferr_t]
labels = ['x', 'y', 'z', 't']
xlab = np.arange(0, t_sink+1, 1)



type_labels = ['Point', '$G1$', '$\\nabla^2 G1$', '$G2$', '$\\nabla^2 G2$', 'Var', 'Var + S/N', 'S/N']
# plt.figure(1, figsize=(16, 12))
# for x in range(5):
#     plt.subplot(511+x)
#     plt.errorbar(xlab, np.resize(rho_mass[x], 32), yerr=np.resize(rho_masserr[x], 32), color='b', ecolor='b', fmt='o', capsize=2)
#     plt.ylabel('$m_{eff}$', fontsize=16)
#     plt.xlabel('$ \Delta t $ (%s)'  % type_labels[x], fontsize=16)
#     plt.xlim(-0.5, 31.5)
#     # plt.yscale('log', nonposy='clip')
#     plt.ylim(0, 2)
#     xp = np.linspace(fitrange[x][0], fitrange[x][1], 100)
#     yp = np.zeros(len(xp))
#     for i in range(len(xp)):
#         yp[i] = masses[x]
#     plt.plot(xp, yp, color='r')

# plt.figure(2, figsize=(16, 12))
# for x in range(3):
#     plt.subplot(311+x)
#     plt.errorbar(xlab, np.resize(rho_mass[x+5], 32), yerr=np.resize(rho_masserr[x+5], 32), color='b', ecolor='b', fmt='o', capsize=2)
#     plt.ylabel('$m_{eff}$', fontsize=16)
#     plt.xlabel('$ \Delta t $ (%s)'  % type_labels[x+5], fontsize=16)
#     plt.xlim(-0.5, 31.5)
#     # plt.yscale('log', nonposy='clip')
#     plt.ylim(0, 2)
#     xp = np.linspace(fitrange[x+5][0], fitrange[x+5][1], 100)
#     yp = np.zeros(len(xp))
#     for i in range(len(xp)):
#         yp[i] = masses[x+5]
#     plt.plot(xp, yp, color='r')

# plt.figure(2, figsize=(15, 12))
# for x in range(3,6):
#     plt.subplot(311+x-3)
#     plt.errorbar(xlab, np.resize(rho_mass[x], 32), yerr=np.resize(rho_masserr[x], 32), color='b', ecolor='b', fmt='^', capsize=2)
#     plt.ylabel('$m_{eff}$ (%s)' % type_labels[x], fontsize=16)
#     plt.xlabel('$ \Delta t $', fontsize=16)
#     plt.xlim(-0.5, 31.5)

plt.figure(1, figsize=(16, 12))
for x in range(5):
    plt.subplot(511+x)
    plt.errorbar(xlab, proton_ff_t[x], yerr=proton_fferr_t[x], color='b', ecolor='b', fmt='o', capsize=2)
    plt.ylabel('FF', fontsize=16)
    plt.xlabel('$ \Delta t $ (%s)'  % type_labels[x], fontsize=16)
    # plt.xlim(-0.5, 31.5)
    # plt.yscale('log', nonposy='clip')
    # plt.ylim(0, 2)
    # xp = np.linspace(fitrange[x+5][0], fitrange[x+5][1], 100)
    # yp = np.zeros(len(xp))
    # for i in range(len(xp)):
    #     yp[i] = masses[x+5]
    # plt.plot(xp, yp, color='r')
plt.subplots_adjust(left=0.07, right=0.97, top=0.96, bottom=0.06, wspace=0.12, hspace=0.43)

plt.figure(2, figsize=(16, 12))
for x in range(5,8):
    plt.subplot(311+x-5)
    plt.errorbar(xlab, proton_ff_t[x], yerr=proton_fferr_t[x], color='b', ecolor='b', fmt='^', capsize=2)
    plt.ylabel('$m_{eff}$ (%s)' % type_labels[x], fontsize=16)
    plt.xlabel('$ \Delta t $', fontsize=16)
    plt.ylim(0.4, 1.3)
    # plt.xlim(-0.5, 31.5)

plt.figure(3, figsize=(16, 8))
plt.errorbar(xlab, proton_ff_t[5], yerr=proton_fferr_t[5], color='b', ecolor='b', fmt='^', capsize=2, label='Variational')
plt.errorbar(xlab+0.1, proton_ff_t[6], yerr=proton_fferr_t[6], color='r', ecolor='r', fmt='v', capsize=2, label='Var + S/N')
plt.errorbar(xlab+0.2, proton_ff_t[7], yerr=proton_fferr_t[7], color='g', ecolor='g', fmt='o', capsize=2, label='S/N')
plt.ylabel('FF', fontsize=20)
plt.xlabel('$ \Delta t $', fontsize=20)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode='expand', borderaxespad=0.)
plt.suptitle('scalar current, $q$=(1,0,0)', fontsize=24)
plt.subplots_adjust(top=0.85)

# plt.figure(3, figsize=(16, 12))
# plt.subplot(211)
# plt.errorbar(xlab, proton_ff_z[0], yerr=proton_fferr_z[0], color='b', ecolor='b', fmt='^', capsize=2)
# plt.ylabel('$FF$', fontsize=16)
# plt.xlabel('$ \\tau$ $(\\gamma_3\\gamma_5)$', fontsize=16)
# # plt.xlim(-0.5, 31.5)
# plt.subplot(212)
# plt.errorbar(xlab, proton_ff_t[0], yerr=proton_fferr_t[0], color='b', ecolor='b', fmt='^', capsize=2)
# plt.ylabel('$FF$', fontsize=16)
# plt.xlabel('$ \\tau$ $(\\gamma_4)$', fontsize=16)

plt.show()

# src = np.matrix([[1],[0],[0],[0]])
# twoptfn_corr = np.zeros((len(rho_2ptfn_snkp_mat), len(rho_2ptfn_snkp_mat[0])), dtype=np.complex128)
# for x in range(len(rho_2ptfn_snkp_mat)):
#     twoptfn_corr[x] = compute_corr(src, src, rho_2ptfn_snkp_mat[x])


# nmin = 8
# nmax= 15

# p = fits.fit_corr(twoptfn_corr, nmin, nmax, [0.0005, 0.5], 'exp')

# plt.figure(1, figsize=(12,8))
# x_lab = np.arange(0, len(rho_mass[0]), 1)
# plt.errorbar(x_lab, rho_mass[0], yerr=rho_masserr[0], color='b', ecolor='b', fmt='^', capsize=2)
# x = np.linspace(nmin-0.5, nmax+0.5, 100)
# y = np.zeros(len(x))
# for i in range(len(x)):
#     y[i] = p[1]
# plt.plot(x, y)
# plt.show()

