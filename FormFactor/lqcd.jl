module lqcd
using ProgressMeter

export read_hadspec_file, read_bar3ptfn_file, make_hermitian!, find_ff_with2ptfn

function read_hadspec_file(had::String, m::String, pf::Tuple{Int64,Int64,Int64}, sources::Array{String,1}, sinks::Array{String,1}, file_prefix::String)
    n_sources::Int64 = length(sources)
    n_sinks::Int64 = length(sinks)
    data::Array{Complex{Float64},4} = zeros(Complex{Float64}, 1, 1, n_sources, n_sinks)
    filename::String = ""
    for x::Int64 in 1:n_sources
        for y::Int64 in 1:n_sinks
            if pf == (0,0,0)
                filename = file_prefix * "hadspec/$(had).D$(m).$(sources[x]).$(sinks[y]).sh_$(x-1)_sh_$(y-1).SS"
            else
                filename = file_prefix * "hadspec/$(had)_px$(pf[1])_py$(pf[2])_pz$(pf[3]).D$(m).$(sources[x]).$(sinks[y]).sh_$(x-1)_sh_$(y-1).SS"
            end
            counter::Int64 = -1
            n_configs::Int64 = 0
            t_size::Int64 = 0
            entries::Array{SubString{String},1} = []
            open(filename) do f::IOStream
                lines::Array{String,1} = readlines(f)
                for line::String in lines
                    entries = split(line)
                    if counter == -1
                        n_configs = parse(Int64, entries[1])
                        t_size = parse(Int64, entries[2])
                        if x == 1 && y == 1
                            data = zeros(Complex{Float64}, n_configs, t_size, n_sources, n_sinks)
                        end
                    else
                        data[Int64(floor(counter/t_size)+1), counter%t_size+1, x, y] = parse(Float64, entries[2]) + parse(Float64, entries[3])im
                    end
                    counter += 1
                end
            end
            if counter != n_configs*t_size
                println("Error: File seems to be corrupted")
            end
        end
    end
    return data
end


function read_bar3ptfn_file(seqsource::String, current::String, i::Int64, g::Int64, q::Tuple{Int64,Int64,Int64}, pf::Tuple{Int64,Int64,Int64}, sources::Array{String,1}, sinks::Array{String,1}, file_prefix::String)
    n_sources::Int64 = length(sources)
    n_sinks::Int64 = length(sinks)
    data::Array{Complex{Float64},4} = zeros(Complex{Float64}, 1, 1, n_sources, n_sinks)
    filename::String = ""
    for x::Int64 in 1:n_sources
        for y::Int64 in 1:n_sinks
            filename = file_prefix * "bar3ptfn/$(current)_cur3ptfn_$(seqsource)_i$(i)_g$(g)_qx$(q[1])_qy$(q[2])_qz$(q[3])_pfx$(pf[1])_pfy$(pf[2])_pfz$(pf[3]).$(sources[x]).$(sinks[y]).sh_$(x-1)_sh_$(y-1).SS"
            counter::Int64 = -1
            n_configs::Int64 = 0
            t_size::Int64 = 0
            entries::Array{SubString{String},1} = []
            open(filename) do f::IOStream
                lines::Array{String,1} = readlines(f)
                for line::String in lines
                    entries = split(line)
                    if counter == -1
                        n_configs = parse(Int64, entries[1])
                        t_size = parse(Int64, entries[2])
                        if x == 1 && y == 1
                            data = zeros(Complex{Float64}, n_configs, t_size, n_sources, n_sinks)
                        end
                    else
                        data[Int64(floor(counter/t_size)+1), counter%t_size+1, x, y] = parse(Float64, entries[2]) + parse(Float64, entries[3])im
                    end
                    counter += 1
                end
            end
            if counter != n_configs*t_size
                println("Error: File seems to be corrupted")
            end
        end
    end
    return data
end


function average(data::Array{Complex{Float64},4})
    n_tot::Int64 = size(data, 1)
    av_data::Array{Complex{Float64},3} = zeros(Complex{Float64}, size(data, 2), size(data, 3), size(data, 4))
    for x::Int64 in 1:n_tot
        av_data += data[x,:,:,:]
    end
    av_data /= Float64(n_tot)
    return av_data
end
function average(data::Array{Float64,3})
    n_tot::Int64 = size(data, 1)
    av_data::Array{Float64,2} = zeros(Float64, size(data, 2), size(data, 3))
    for x::Int64 in 1:n_tot
        av_data += data[x,:,:]
    end
    av_data /= Float64(n_tot)
    return av_data
end


function make_hermitian!(data::Array{Complex{Float64},4})
    n_configs::Int64 = size(data,1)
    t_size::Int64 = size(data,2)
    n_sources::Int64 = size(data,3)
    n_sinks::Int64 = size(data,4)
    temp_re::Float64 = 0.
    temp_im::Float64 = 0.
    for n::Int64 in 1:n_configs
        if n == n_configs
            println("\rMaking matrices hermitian... \033[1;32mdone\033[0m")
        else
            print("\rMaking matrices hermitian... ", floor(n/n_configs*100.),"%")
        end
        for t::Int64 in 1:t_size
            for x::Int64 in 1:n_sources
                for y::Int64 in 1:n_sinks
                    if x == y
                        data[n, t, x, y] = real(data[n, t, x, y])
                    else
                        temp_re = (real(data[n, t, x, y]) + real(data[n, t, y, x]))/2.0
                        temp_im = (imag(data[n, t, x, y]) - imag(data[n, t, y, x]))/2.0
                        data[n, t, x, y] = temp_re + temp_im*1.0im
                        data[n, t, y, x] = temp_re - temp_im*1.0im
                    end
                end
            end
        end
    end
    return data
end


function order!(evals::Array{Complex{Float64},1}, evecs::Array{Complex{Float64},2})
    temp_evals::Array{Complex{Float64},1} = copy(evals)
    temp_evecs::Array{Complex{Float64},2} = copy(evecs)
    ord::Array{Int64,1} = sortperm(abs(evals))
    for i::Int64 in 1:length(ord)
        evals[i] = temp_evals[ord[i]]
        evecs[:,i] = temp_evecs[:,ord[i]]
    end
    return evals, evecs
end


function normalize_evecs!(evecs::Array{Complex{Float64},2})
    for i::Int64 in 1:length(evecs[1,:])
        evecs[:,i] = evecs[:,i]/norm(evecs[:,i])
    end
    return evecs
end


function find_eigsys(data::Array{Complex{Float64},3}, t_gen_ev::Int64)
    t_size::Int64 = size(data,1)
    n_sources::Int64 = size(data,2)
    n_sinks::Int64 = size(data,3)
    ## IMPORTANT
    # since Julia starts indices with 1 then we must
    # add 1 to t_gen_ev
    if n_sources != n_sinks
        println("Error: Matrix is not square")
    end
    evals::Array{Complex{Float64},2} = zeros(Complex{Float64}, t_size, n_sources)
    evecs::Array{Complex{Float64},3} = zeros(Complex{Float64}, t_size, n_sources, n_sources)
    for t::Int64 in 1:t_size
        evals[t,:], evecs[t,:,:] = eig(data[t,:,:], data[t_gen_ev+1,:,:])
        order!(evals[t,:], evecs[t,:,:])
        normalize_evecs!(evecs[t,:,:])
    end
    return evals, evecs
end


function find_sn_sink(source::Array{Complex{Float64},1}, mat::Array{Complex{Float64},4}, t_sn::Int64)
    n_configs::Int64 = size(mat, 1)
    n_sources::Int64 = size(mat, 3)
    n_sinks::Int64 = size(mat, 4)
    if n_sources != n_sinks
        print("Error: Matrix is not square")
    end
    ## IMPORTANT
    # since Julia starts indices with 1 then we must
    # add 1 to t_sn
    av_mat::Array{Complex{Float64},3} = average(mat)
    sink::Array{Complex{Float64},1} = zeros(Complex{Float64}, n_sources)
    sigma2::Array{Complex{Float64},2} = zeros(Complex{Float64}, n_sources, n_sources)
    for n::Int64 in 1:n_configs
        sigma2 += mat[n,t_sn+1,:,:] * reshape(source, n_sources, 1) * reshape(source, n_sources, 1)' * mat[n,t_sn+1,:,:]'
    end
    sigma2 /= Float64(n_configs)
    inv_sigma2::Array{Complex{Float64},2} = inv(sigma2)
    a = (reshape(source, n_sources, 1)' * av_mat[t_sn+1,:,:]' * inv_sigma2 * inv_sigma2 * av_mat[t_sn+1,:,:] * reshape(source, n_sources, 1))[1,1]
    a = 1.0/sqrt(a)
    sink = (a * inv_sigma2 * av_mat[t_sn+1,:,:] * reshape(source, n_sources, 1))[:,1]
    sink /= sqrt((reshape(sink, n_sources, 1)' * reshape(sink, n_sources, 1))[1,1])
    return sink
end


function find_opt_sn(mat::Array{Complex{Float64},4}, t_sn::Int64)
    n_configs::Int64 = size(mat, 1)
    n_sources::Int64 = size(mat, 3)
    n_sinks::Int64 = size(mat, 4)
    if n_sources != n_sinks
        print("Error: Matrix is not square")
    end
    av_mat::Array{Complex{Float64},3} = average(mat)
    n_iter::Int64 = 20
    n_tries::Int64 = 5
    converged::Bool = false
    source::Array{Complex{Float64},2} = zeros(Complex{Float64}, n_sources, 1)
    sink::Array{Complex{Float64},2} = zeros(Complex{Float64}, n_sources, 1)
    ## IMPORTANT
    # since Julia starts indices with 1 then we must
    # add 1 to t_sn
    for tries::Int64 in 1:n_tries
        source = reshape([1.+0.0im, 0.0im,0.0im,0.0im,0.0im],5,1) #reshape(randn(n_sources) + randn(n_sources)im, n_sources, 1)
        source /= norm(source)
        sink = copy(source)
        old_source::Array{Complex{Float64},2} = copy(source)
        old_sink::Array{Complex{Float64},2} = copy(sink)
        for n::Int64 in 1:n_iter
            sigma2_source::Array{Complex{Float64},2} = zeros(Complex{Float64}, n_sources, n_sources)
            sigma2_sink::Array{Complex{Float64},2} = zeros(Complex{Float64}, n_sources, n_sources)
            for x::Int64 in 1:n_configs
                sigma2_source += mat[x,t_sn+1,:,:] * reshape(source, n_sources, 1) * reshape(source, n_sources, 1)' * mat[x,t_sn+1,:,:]'
                sigma2_sink += mat[x,t_sn+1,:,:]' * reshape(sink, n_sources, 1) * reshape(sink, n_sources, 1)' * mat[x,t_sn+1,:,:]
            end
            sigma2_source /= Float64(n_configs)
            sigma2_sink /= Float64(n_configs)
            inv_sigma2_source::Array{Complex{Float64},2} = inv(sigma2_source)
            inv_sigma2_sink::Array{Complex{Float64},2} = inv(sigma2_sink)
            a_source = (reshape(source, n_sources, 1)' * av_mat[t_sn+1,:,:]' * inv_sigma2_source * inv_sigma2_source * av_mat[t_sn+1,:,:] * reshape(source, n_sources, 1))[1,1]
            a_sink = (reshape(sink, n_sources, 1)' * av_mat[t_sn+1,:,:] * inv_sigma2_sink * inv_sigma2_sink * av_mat[t_sn+1,:,:]' * reshape(sink, n_sources, 1))[1,1]
            old_source = copy(source)
            old_sink = copy(sink)
            source = a_sink * inv_sigma2_sink * av_mat[t_sn+1,:,:]' * reshape(old_sink, n_sources, 1)
            sink = a_source * inv_sigma2_source * av_mat[t_sn+1,:,:] * reshape(old_source, n_sources, 1)
            source /= norm(source)
            sink /= norm(sink)
        end
        if 1 - real((reshape(old_source, n_sources, 1)' * reshape(old_source, n_sources, 1))[1,1]) < 0.01
            converged = true
            break
        end
    end
    if !converged
        println("S/N process did not converge")
        return zeros(Complex{Float64}, n_sources, 1), zeros(Complex{Float64}, n_sources, 1)
    end
    #println(source)
    return reshape(source, n_sources), reshape(sink, n_sources)
end


function compute_corr(source::Array{Complex{Float64},1}, sink::Array{Complex{Float64},1}, mat::Array{Complex{Float64},3})
    t_size::Int64 = size(mat, 1)
    n_sources::Int64 = size(mat, 2)
    n_sinks::Int64 = size(mat, 3)
    corr::Array{Complex{Float64},1} = zeros(Complex{Float64}, t_size)
    for x::Int64 in 1:t_size
        corr[x] = (reshape(sink, n_sinks, 1)' * mat[x,:,:] * reshape(source, n_sources, 1))[1,1]
    end
    return corr
end


function compute_ff_complex(threeptfn::Array{Complex{Float64},1}, twoptfn_src::Array{Complex{Float64},1}, twoptfn_snk::Array{Complex{Float64},1}, t_sink::Int64)
    ff::Array{Float64,1} = zeros(Float64, t_sink+1)
    r::Complex{Float64} = 0.
    ## IMPORTANT
    # since Julia starts indices with 1 then we must
    # add 1 to t_sink
    fixed_t_sink::Int64 = t_sink+1
    for tau::Int64 in 1:t_sink+1
        r = twoptfn_snk[fixed_t_sink]/twoptfn_src[fixed_t_sink]
        r *= twoptfn_snk[tau]/twoptfn_src[tau]
        r *= twoptfn_src[fixed_t_sink-tau+1]/twoptfn_snk[fixed_t_sink-tau+1]
        r = sqrt(r)
        r *= threeptfn[tau]/twoptfn_snk[fixed_t_sink]
        ff[tau] = abs(r)
    end
    return ff
end


function find_ff_with2ptfn(threeptfn::Array{Complex{Float64},4}, twoptfn_src::Array{Complex{Float64},4}, twoptfn_snk::Array{Complex{Float64},4}, t_sink::Int64, t_gen_ev::Int64, t_var::Int64, t_sn::Int64; n_boot::Int64 = 100)
    p = Progress(n_boot, 0.5, "Performing bootstrapping... ", 50)
    n_configs::Int64 = size(threeptfn, 1)
    t_size::Int64 = size(threeptfn, 2)
    n_sources::Int64 = size(threeptfn, 3)
    n_sinks::Int64 = size(threeptfn, 4)
    n_failed::Int64 = 0
    temp_rnd_conf::Array{Int64,1} = []
    all_ff::Array{Float64,3} = zeros(Float64, n_boot, 8, t_sink+1)
    if n_sources != n_sinks
        println("Error: Matrix is not square")
    end
    for n::Int64 in 1:n_boot
        #if n == n_boot
        #    println("\rPerforming bootstrapping for charge/form factor... \033[1;32mdone\033[0m")
        #else
        #    print("\rPerforming bootstrapping for charge/form factor... ", Int64(floor(n/n_boot*100)), "%")
        #end
        temp_3ptfn::Array{Complex{Float64},4} = zeros(Complex{Float64}, n_configs, t_size, n_sources, n_sources)
        temp_2ptfn_src::Array{Complex{Float64},4} = zeros(Complex{Float64}, n_configs, t_size, n_sources, n_sources)
        temp_2ptfn_snk::Array{Complex{Float64},4} = zeros(Complex{Float64}, n_configs, t_size, n_sources, n_sources)
        temp_rnd_conf = rand(1:n_configs, n_configs)
        for i::Int64 in 1:n_configs
            temp_3ptfn[i,:,:,:] = threeptfn[temp_rnd_conf[i],:,:,:]
            temp_2ptfn_src[i,:,:,:] = twoptfn_src[temp_rnd_conf[i],:,:,:]
            temp_2ptfn_snk[i,:,:,:] = twoptfn_snk[temp_rnd_conf[i],:,:,:]
        end
        av_3ptfn::Array{Complex{Float64},3} = average(temp_3ptfn)
        av_2ptfn_src::Array{Complex{Float64},3} = average(temp_2ptfn_src)
        av_2ptfn_snk::Array{Complex{Float64},3} = average(temp_2ptfn_snk)

        ## Smeared src and snk ##
        id_matrix::Array{Complex{Float64},2} = eye(Complex{Float64}, n_sources)
        for s::Int64 in 1:n_sources
            smeared_3ptfn_corr::Array{Complex{Float64},1} = compute_corr(id_matrix[:,s], id_matrix[:,s], av_3ptfn)
            smeared_2ptfn_src_corr::Array{Complex{Float64},1} = compute_corr(id_matrix[:,s], id_matrix[:,s], av_2ptfn_src)
            smeared_2ptfn_snk_corr::Array{Complex{Float64},1} = compute_corr(id_matrix[:,s], id_matrix[:,s], av_2ptfn_snk)
            all_ff[n, s, :] = compute_ff_complex(smeared_3ptfn_corr, smeared_2ptfn_src_corr, smeared_2ptfn_snk_corr, t_sink)
        end

        ## Variational src and snk ##
        evals::Array{Complex{Float64},2}, evecs::Array{Complex{Float64},3} = find_eigsys(av_2ptfn_snk, t_gen_ev)
        var_src::Array{Complex{Float64},1} = evecs[t_var+1,:,1]
        var_3ptfn_corr::Array{Complex{Float64},1} = compute_corr(var_src, var_src, av_3ptfn)
        var_2ptfn_src_corr::Array{Complex{Float64},1} = compute_corr(var_src, var_src, av_2ptfn_src)
        var_2ptfn_snk_corr::Array{Complex{Float64},1} = compute_corr(var_src, var_src, av_2ptfn_snk)
        all_ff[n, 6, :] = compute_ff_complex(var_3ptfn_corr, var_2ptfn_src_corr, var_2ptfn_snk_corr, t_sink)

        ## Var source and S/N sink ##
        varsn_snk::Array{Complex{Float64},1} = find_sn_sink(var_src, temp_2ptfn_snk, t_sn)
        varsn_3ptfn_corr::Array{Complex{Float64},1} = compute_corr(var_src, varsn_snk, av_3ptfn)
        varsn_2ptfn_src_corr::Array{Complex{Float64},1} = compute_corr(var_src, varsn_snk, av_2ptfn_src)
        varsn_2ptfn_snk_corr::Array{Complex{Float64},1} = compute_corr(var_src, varsn_snk, av_2ptfn_snk)
        all_ff[n, 7, :] = compute_ff_complex(varsn_3ptfn_corr, varsn_2ptfn_src_corr, varsn_2ptfn_snk_corr, t_sink)

        ## S/N source and sink ##
        sn_src::Array{Complex{Float64},1}, sn_snk::Array{Complex{Float64},1} = find_opt_sn(temp_2ptfn_snk, t_sn)
        if sn_src == zeros(Complex{Float64}, n_sources)
            n_failed += 1
        else
            sn_2ptfn_src_corr::Array{Complex{Float64},1} = compute_corr(sn_src, sn_snk, av_2ptfn_src)
            sn_2ptfn_snk_corr::Array{Complex{Float64},1} = compute_corr(sn_src, sn_snk, av_2ptfn_snk)
            sn_3ptfn_corr::Array{Complex{Float64},1} = compute_corr(sn_src, sn_snk, av_3ptfn)
            all_ff[n, 8, :] = compute_ff_complex(sn_3ptfn_corr, sn_2ptfn_src_corr, sn_2ptfn_snk_corr, t_sink)
        end

        ## update progress bar ##
        update!(p, n)
    end
    av_ff::Array{Float64,2} = average(all_ff)
    if n_failed != n_boot
        av_ff[8,:] *= Float64(n_boot)/Float64(n_boot-n_failed)
    end
    err_ff::Array{Float64,2} = zeros(Float64, 8, t_sink+1)
    for n::Int64 in 1:n_boot
        err_ff += (all_ff[n, :, :]-av_ff).^2
    end
    err_ff = sqrt(err_ff/Float64(n_boot))
    if n_failed != n_boot
        err_ff[8,:] *= sqrt(Float64(n_boot)/Float64(n_boot-n_failed))
    end
    return av_ff, err_ff, all_ff
end


end
