module lqcdfits
using Optim, lqcd, ProgressMeter

export exp_func, cosh_func, sinh_func, const_func, fit_corr, find_2ptcorr_wmat, scan_range_2ptfn, find_mass_fit

function exp_func(t::Int64, Nt::Int64, par::Array{Float64,1})
    return par[1]*exp(-par[2]*t)
end


function cosh_func(t::Int64, Nt::Int64, par::Array{Float64,1})
    return par[1]*cosh((t-Nt/2.)*par[2])
end


function sinh_func(t::Int64, Nt::Int64, par::Array{Float64,1})
    return par[1]*sinh((t-Nt/2.)*par[2])
end


function const_func(t::Int64, Nt::Int64, par::Float64)
    return par
end


function weight_mat(corr::Array{Float64,2})
    n_configs::Int64 = size(corr, 1)
    Nt::Int64 = size(corr, 2)
    av_corr::Array{Float64,1} = average(corr)
    cov_mat::Array{Float64,2} = zeros(Float64, Nt, Nt)
    for x::Int64 in 1:Nt
        for y::Int64 in 1:Nt
            for config::Int64 in 1:n_configs
                cov_mat[x,y] += (corr[config,x] - av_corr[x])*(corr[config,y] - av_corr[y])
            end
        end
    end
    cov_mat /= Float64(n_configs*(n_configs-1))
    w_mat::Array{Float64,2} = inv(cov_mat)
    return w_mat
end


function chi2(param, av_corr::Array{Float64,1}, w_mat::Array{Float64,2}, fit_func, Nt::Int64, nmin::Int64, nmax::Int64)
    x2::Float64 = 0.
    for x::Int64 in nmin:nmax
        for y::Int64 in nmin:nmax
            x2 += (av_corr[x] - fit_func(x,Nt,param)) * w_mat[x,y] * (av_corr[y] - fit_func(y,Nt,param))
        end
    end
    x2 /= nmax-nmin+1-size(param,1)
    return x2
end


function fit_corr(corr::Array{Float64,2}, nmin::Int64, nmax::Int64, par_guess::Array{Float64,1}, func::String="const")
    if func == "exp"
        fit_func = exp_func
    elseif func == "cosh"
        fit_func = cosh_func
    elseif func == "sinh"
        fit_func = sinh_func
    elseif func == "const"
        fit_func = const_func
    else
        println("Error: Function type not defined, using const")
        fit_func = const_func
    end
    av_corr::Array{Float64,1} = average(corr)
    Nt::Int64 = size(av_corr, 1)
    w_mat::Array{Float64,2} = weight_mat(corr)
    function tmp_function(par)
        return chi2(par, av_corr, w_mat, fit_func, Nt, nmin, nmax)
    end
    if func == "const"
        res = optimize(tmp_function, par_guess[1], par_guess[2])
    else
        res = optimize(tmp_function, par_guess)
    end
    if Optim.converged(res)
        #println("\033[1;32mMinimization succeeded\033[0m")
        #println(Optim.minimum(res))
    else
        print("\033[1;31mMinimization failed\033[0m")
    end
    return Optim.minimizer(res), Optim.minimum(res)
end
function fit_corr(av_corr::Array{Float64,1}, w_mat::Array{Float64,2}, nmin::Int64, nmax::Int64, par_guess::Array{Float64,1}, func::String="const")
    if func == "exp"
        fit_func = exp_func
    elseif func == "cosh"
        fit_func = cosh_func
    elseif func == "sinh"
        fit_func = sinh_func
    elseif func == "const"
        fit_func = const_func
    else
        println("Error: Function type not defined, using const")
        fit_func = const_func
    end
    Nt::Int64 = size(av_corr, 1)
    function tmp_function(par)
        return chi2(par, av_corr, w_mat, fit_func, Nt, nmin, nmax)
    end
    if func == "const"
        res = optimize(tmp_function, par_guess[1], par_guess[2])
    else
        res = optimize(tmp_function, par_guess)
    end
    if Optim.converged(res)
        #println("\033[1;32mMinimization succeeded\033[0m")
        #println(Optim.minimizer(res))
    else
        print("\033[1;31mMinimization failed\033[0m")
    end
    return Optim.minimizer(res), Optim.minimum(res)
end


function find_2ptcorr_wmat(twoptfn::Array{Complex{Float64},4}, t_gen_ev::Int64, t_var::Int64, t_sn::Int64)
    n_configs::Int64 = size(threeptfn, 1)
    t_size::Int64 = size(threeptfn, 2)
    n_sources::Int64 = size(threeptfn, 3)
    n_sinks::Int64 = size(threeptfn, 4)
    w_mats::Array{Float64,3} = zeros(Float64, 8, n_sources, n_sources)
    if n_sources != n_sinks
        println("Error: Matrix is not square")
    end
    av_2ptfn::Array{Complex{Float64},3} = average(temp_2ptfn)

    ## Smeared src and snk ##
    id_matrix::Array{Complex{Float64},2} = eye(Complex{Float64}, n_sources)
    for s::Int64 in 1:n_sources
        w_mats[s, :, :] = weight_mat(real(compute_corrs(id_matrix[:,s], id_matrix[:,s], twoptfn)))
    end

    ## Variational src and snk ##
    evals::Array{Complex{Float64},2}, evecs::Array{Complex{Float64},3} = find_eigsys(av_2ptfn, t_gen_ev)
    var_src::Array{Complex{Float64},1} = evecs[t_var+1,:,1]
    w_mats[6, :, :] = weight_mat(real(compute_corrs(var_src, var_src, twoptfn)))

    ## Var source and S/N sink ##
    varsn_snk::Array{Complex{Float64},1} = find_sn_sink(var_src, temp_2ptfn_snk, t_sn)
    w_mats[7, :, :] = weight_mat(real(compute_corrs(var_src, varsn_snk, twoptfn)))

    ## S/N source and sink ##
    sn_src::Array{Complex{Float64},1}, sn_snk::Array{Complex{Float64},1} = find_opt_sn(temp_2ptfn_snk, t_sn)
    if sn_src != zeros(Complex{Float64}, n_sources)
        w_mats[8, :, :] = weight_mat(real(compute_corrs(sn_src, sn_snk, twoptfn)))
    end

    return w_mats
end


function scan_range_2ptfn(twoptfn::Array{Complex{Float64},4}, t_gen_ev::Int64, t_var::Int64, t_sn::Int64, par_guess::Array{Float64,1}; nmin::Int64=10, nmax::Int64=17, func::String="exp")
    n_configs::Int64 = size(twoptfn, 1)
    t_size::Int64 = size(twoptfn, 2)
    n_sources::Int64 = size(twoptfn, 3)
    n_sinks::Int64 = size(twoptfn, 4)
    w_mats::Array{Float64,3} = zeros(Float64, 8, t_size, t_size)
    av_corrs::Array{Float64,2} = zeros(Float64, 8, t_size)
    if n_sources != n_sinks
        println("Error: Matrix is not square")
    end
    av_2ptfn::Array{Complex{Float64},3} = average(twoptfn)

    ## Smeared src and snk ##
    id_matrix::Array{Complex{Float64},2} = eye(Complex{Float64}, n_sources)
    for s::Int64 in 1:n_sources
        w_mats[s, :, :] = weight_mat(real(compute_corrs(id_matrix[:,s], id_matrix[:,s], twoptfn)))
        av_corrs[s,:] = real(compute_corr(id_matrix[:,s], id_matrix[:,s], av_2ptfn))
    end

    ## Variational src and snk ##
    evals::Array{Complex{Float64},2}, evecs::Array{Complex{Float64},3} = find_eigsys(av_2ptfn, t_gen_ev)
    var_src::Array{Complex{Float64},1} = evecs[t_var+1,:,1]
    w_mats[6, :, :] = weight_mat(real(compute_corrs(var_src, var_src, twoptfn)))
    av_corrs[6,:] = real(compute_corr(var_src, var_src, av_2ptfn))

    ## Var source and S/N sink ##
    varsn_snk::Array{Complex{Float64},1} = find_sn_sink(var_src, twoptfn, t_sn)
    w_mats[7, :, :] = weight_mat(real(compute_corrs(var_src, varsn_snk, twoptfn)))
    av_corrs[7,:] = real(compute_corr(var_src, varsn_snk, av_2ptfn))

    ## S/N source and sink ##
    sn_src::Array{Complex{Float64},1}, sn_snk::Array{Complex{Float64},1} = find_opt_sn(twoptfn, t_sn)
    if sn_src != zeros(Complex{Float64}, n_sources)
        w_mats[8, :, :] = weight_mat(real(compute_corrs(sn_src, sn_snk, twoptfn)))
        av_corrs[8,:] = real(compute_corr(sn_src, sn_snk, av_2ptfn))
    end

    limits::Array{Int64,2} = zeros(Int64,8,2)
    mass_sys_err::Array{Float64,1} = zeros(Float64,8)
    mass_sys::Array{Float64,1} = zeros(Float64,8)
    p = Progress(8, 1, "Scanning fit range...         ", 50)
    for i::Int64 in 1:8
        n_low::Array{Int64,1} = zeros(Int64,0)
        n_high::Array{Int64,1} = zeros(Int64,0)
        chi_low::Array{Float64,1} = zeros(Float64,0)
        chi_high::Array{Float64,1} = zeros(Float64,0)
        for x::Int64 in 4:(nmax-1)
            chi2::Float64 = fit_corr(av_corrs[i,:], w_mats[i,:,:], nmax-x, nmax, par_guess, func)[2]
            push!(n_low, nmax-x)
            push!(chi_low, chi2)
            if chi2 > 10
                break
            end
        end
        for x::Int64 in 4:(t_size-nmin)
            chi2::Float64 = fit_corr(av_corrs[i,:], w_mats[i,:,:], nmin, nmin+x, par_guess, func)[2]
            push!(n_high, nmin+x)
            push!(chi_high, chi2)
            if chi2 > 10
                break
            end
        end

        n_min::Int64 = 0
        n_max::Int64 = 0
        for x::Int64 in 1:length(n_low)
            if chi_low[x] < 3
                n_min = n_low[x]
            end
        end
        for x::Int64 in 1:length(n_high)
            if chi_high[x] < 3 && n_high[x] < t_size/2
                n_max = n_high[x]
            end
        end

        limits[i,1] = n_min
        limits[i,2] = n_max
        all_mass::Array{Float64,1} = zeros(Float64,0)

        for x::Int64 in -1:1
            for y::Int64 in -1:1
                mass::Float64 = fit_corr(av_corrs[i,:], w_mats[i,:,:], n_min+x, n_max+y, par_guess, func)[1][2]
                push!(all_mass, mass)
            end
        end

        av_mass::Float64 = mean(all_mass)
        mass_sys[i] = av_mass
        std_mass::Float64 = 0.
        for x::Int64 in 1:length(all_mass)
            std_mass += (all_mass[x] - av_mass)^2
        end
        std_mass = sqrt(std_mass)/9.
        mass_sys_err[i] = std_mass

        update!(p, i)
    end
    return limits, mass_sys, mass_sys_err
end


function find_mass_fit(twoptfn::Array{Complex{Float64},4}, t_gen_ev::Int64, t_var::Int64, t_sn::Int64, limits::Array{Int64,2}, par_guess::Array{Float64,1}; func::String="exp", n_boot::Int64 = 100)
    p = Progress(n_boot, 0.5, "Performing bootstrapping...   ", 50)
    n_configs::Int64 = size(twoptfn, 1)
    t_size::Int64 = size(twoptfn, 2)
    n_sources::Int64 = size(twoptfn, 3)
    n_sinks::Int64 = size(twoptfn, 4)
    n_failed::Int64 = 0
    temp_rnd_conf::Array{Int64,1} = []
    all_masses::Array{Float64,2} = zeros(Float64, n_boot, 8)
    if n_sources != n_sinks
        println("Error: Matrix is not square")
    end
    for n::Int64 in 1:n_boot
        temp_2ptfn::Array{Complex{Float64},4} = zeros(Complex{Float64}, n_configs, t_size, n_sources, n_sources)
        temp_rnd_conf = rand(1:n_configs, n_configs)
        for i::Int64 in 1:n_configs
            temp_2ptfn[i,:,:,:] = twoptfn[temp_rnd_conf[i],:,:,:]
        end
        av_2ptfn::Array{Complex{Float64},3} = average(temp_2ptfn)

        ## Smeared src and snk ##
        id_matrix::Array{Complex{Float64},2} = eye(Complex{Float64}, n_sources)
        for s::Int64 in 1:n_sources
            smeared_2ptfn_corr::Array{Float64,2} = real(compute_corrs(id_matrix[:,s], id_matrix[:,s], temp_2ptfn))
            all_masses[n, s] = fit_corr(smeared_2ptfn_corr, limits[s,1], limits[s,2], par_guess, func)[1][2]
        end

        ## Variational src and snk ##
        evals::Array{Complex{Float64},2}, evecs::Array{Complex{Float64},3} = find_eigsys(av_2ptfn, t_gen_ev)
        var_src::Array{Complex{Float64},1} = evecs[t_var+1,:,1]
        var_2ptfn_corr::Array{Float64,2} = real(compute_corrs(var_src, var_src, temp_2ptfn))
        all_masses[n, 6] = fit_corr(var_2ptfn_corr, limits[6,1], limits[6,2], par_guess, func)[1][2]

        ## Var source and S/N sink ##
        varsn_snk::Array{Complex{Float64},1} = find_sn_sink(var_src, temp_2ptfn, t_sn)
        varsn_2ptfn_corr::Array{Float64,2} = real(compute_corrs(var_src, varsn_snk, temp_2ptfn))
        all_masses[n, 7] = fit_corr(varsn_2ptfn_corr, limits[7,1], limits[7,2], par_guess, func)[1][2]

        ## S/N source and sink ##
        sn_src::Array{Complex{Float64},1}, sn_snk::Array{Complex{Float64},1} = find_opt_sn(temp_2ptfn, t_sn)
        if sn_src == zeros(Complex{Float64}, n_sources)
            n_failed += 1
        else
            sn_2ptfn_corr::Array{Float64,2} = real(compute_corrs(sn_src, sn_snk, temp_2ptfn))
            all_masses[n, 8] = fit_corr(sn_2ptfn_corr, limits[8,1], limits[8,2], par_guess, func)[1][2]
        end

        ## update progress bar ##
        update!(p, n)
    end
    av_masses::Array{Float64,1} = average(all_masses)
    if n_failed != n_boot
        av_masses[8] *= Float64(n_boot)/Float64(n_boot-n_failed)
    end
    err_masses::Array{Float64,1} = zeros(Float64, 8)
    for n::Int64 in 1:n_boot
        err_masses += (all_masses[n, :]-av_masses).^2
    end
    err_masses = sqrt(err_masses/Float64(n_boot))
    if n_failed != n_boot
        err_masses[8] *= sqrt(Float64(n_boot)/Float64(n_boot-n_failed))
    end
    return av_masses, err_masses
end


end
