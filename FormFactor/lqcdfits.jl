module lqcdfits
using Optim

export fit_corr

function average(data::Array{Float64,2})
    n_tot::Int64 = size(data, 1)
    av_data::Array{Float64,1} = zeros(Float64, size(data, 2))
    for x::Int64 in 1:n_tot
        av_data += data[x,:]
    end
    av_data /= Float64(n_tot)
    return av_data
end


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
        println("\033[1;32mMinimization succeeded\033[0m")
        println(Optim.minimum(res))
    else
        print("\033[1;31mMinimization failed\033[0m")
    end
    return Optim.minimizer(res), Optim.minimum(res)
end




end
