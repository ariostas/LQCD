module lqcd

const var_type = Complex{Float64}
const storage_type = UInt16
const distance_type = Float32
const index_type = UInt32

const sigma_0 = [one(var_type) zero(var_type);
                 zero(var_type) one(var_type)]
const sigma_1 = [zero(var_type) one(var_type);
                 one(var_type) zero(var_type)]
const sigma_2 = [zero(var_type) -one(var_type)*im;
                 one(var_type)*im zero(var_type)]
const sigma_3 = [one(var_type) zero(var_type);
                 zero(var_type) -one(var_type)]
const sigma_n = [sigma_1,sigma_2,sigma_3]

# Gell-Mann matrices
const T_0 = [one(var_type) zero(var_type) zero(var_type);
             zero(var_type) one(var_type) zero(var_type);
             zero(var_type) zero(var_type) one(var_type)]
const T_1 = [zero(var_type) one(var_type) zero(var_type);
             one(var_type) zero(var_type) zero(var_type);
             zero(var_type) zero(var_type) zero(var_type)]
const T_2 = [zero(var_type) -one(var_type)*im zero(var_type);
             one(var_type)*im zero(var_type) zero(var_type);
             zero(var_type) zero(var_type) zero(var_type)]
const T_3 = [one(var_type) zero(var_type) zero(var_type);
             zero(var_type) -one(var_type) zero(var_type);
             zero(var_type) zero(var_type) zero(var_type)]
const T_4 = [zero(var_type) zero(var_type) one(var_type);
             zero(var_type) zero(var_type) zero(var_type);
             one(var_type) zero(var_type) zero(var_type)]
const T_5 = [zero(var_type) zero(var_type) -one(var_type)*im;
             zero(var_type) zero(var_type) zero(var_type);
             one(var_type)*im zero(var_type) zero(var_type)]
const T_6 = [zero(var_type) zero(var_type) zero(var_type);
             zero(var_type) zero(var_type) one(var_type);
             zero(var_type) one(var_type) zero(var_type)]
const T_7 = [zero(var_type) zero(var_type) zero(var_type);
             zero(var_type) zero(var_type) -one(var_type)*im;
             zero(var_type) one(var_type)*im zero(var_type)]
const T_8 = var_type(1./sqrt(3.))*[one(var_type) zero(var_type) zero(var_type);
                                  zero(var_type) one(var_type) zero(var_type);
                                  zero(var_type) zero(var_type) -var_type(2)]
const T_n = [T_1,T_2,T_3,T_4,T_5,T_6,T_7,T_8]

const μ₁ = 1./2.*(-1.+sqrt(5))
const μ₂ = 1./2.*(-1.-sqrt(5))
const ω = exp(2π*im/3.)
const A = Array{var_type}(
          [0. 1. 0.;
           0. 0. 1.;
           1. 0. 0.])
const B = Array{var_type}(
          [1. 0. 0.;
           0. -1. 0.;
           0. 0. -1.])
const C = Array{var_type}(
          1./2.*[-1. μ₂ μ₁;
                 μ₂ μ₁ -1.;
                 μ₁ -1. μ₂])
const D = Array{var_type}(
          [-1. 0. 0.;
           0. 0. -ω;
           0. -ω^2 0.])
const gen = [eye(var_type,3,3),A,B,C,D]

const id = eye(var_type,3)
const ze = zeros(id)

export A,B,C,D,gen,T_1,T_2,T_3,T_4,T_5,T_6,T_7,T_8,T_n,id,ze,
       fix_unitarity, fix_unitarity!, find_gp_matrices, load_gp_matrices,
       distance, all_bit_combinations, random_matrices, load_random_matrices,
       fixed_points, load_fixed_points, set_up_random_matrices,
       subdiv_gp_matrices, multilateration, distances, multilateration_simd,
       find_closest_matrix_multi, find_closest_matrix_direct


function set_up_random_matrices()
    eps = [0.007, 0.017, 0.027, 0.04]
    for i in eps
        gen_mat, d, m = random_matrices(i)
        fixed_points(gen_mat,i)
    end
    return
end


function all_bit_combinations(n_bits)
    combs = zeros(storage_type,2^n_bits)
    combs[:] = [0:2^n_bits-1;]
    return combs
end


function random_matrices(epsilon, n_bits=16)
    rand_m = Array{var_type}(3,3,n_bits)
    all_combs = all_bit_combinations(n_bits)
    for i in 1:n_bits
        r = rand(8)-0.5
        while(norm(r) > 1)
            r = 2.0*rand(8)-1.
        end
        r *= epsilon
        rand_m[:,:,i] = expm(im*sum(T_n.*r))
        #rand_m[:,:,i] = expm(im*epsilon*sum((rand(8)-0.5).*T_n))
    end
    tmp_dist = map(x->distance(id,x),
                        [prod([rand_m[:,:,j]^(div(all_combs[i],2^(j-1))%2)
                                    for j in 1:n_bits]) for i in 1:2^n_bits])

    order = sortperm(tmp_dist)
    dist = Array{distance_type}(2^n_bits)
    matrices = Array{storage_type}(2^n_bits)
    for i in 1:2^n_bits
        dist[i] = tmp_dist[order[i]]
        matrices[i] = all_combs[order[i]]
    end
    write("Matrix_data/random_matrices_eps_"*string(n_bits)*"bits_"
                            *string(epsilon)*"eps_genmatrices", rand_m)
    write("Matrix_data/random_matrices_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_distances", dist)
    write("Matrix_data/random_matrices_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_matrices", matrices)
    return rand_m, dist, matrices
end


function load_random_matrices(epsilon, n_bits=16)
    dist = Array{distance_type}(2^n_bits)
    matrices = Array{storage_type}(2^n_bits)
    gen_matrices = Array{var_type}(3,3,n_bits)
    read!("Matrix_data/random_matrices_eps_"*string(n_bits)*"bits_"
                            *string(epsilon)*"eps_genmatrices", gen_matrices)
    read!("Matrix_data/random_matrices_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_distances", dist)
    read!("Matrix_data/random_matrices_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_matrices", matrices)
    return gen_matrices, dist, matrices
end


function fixed_points(rand_m, epsilon, epsilon_fixed=0.75)
    n_bits = size(rand_m,3)
    fixed_pts = Array{var_type}(3,3,9)
    for i in 1:8
        fixed_pts[:,:,i] = expm(im*epsilon_fixed*T_n[i])
    end
    fixed_pts[:,:,9] = id
    all_combs = all_bit_combinations(n_bits)
    tmp_dist = Array{distance_type}(2^n_bits,9)
    dist = Array{distance_type}(2^n_bits,9)
    matrices = Array{storage_type}(2^n_bits,9)
    for i in 1:9
        tmp_dist[:,i] = map(x->distance(fixed_pts[:,:,i],x),
                            [prod([rand_m[:,:,j]^(div(all_combs[i],2^(j-1))%2)
                                for j in 1:n_bits]) for i in 1:2^n_bits])

        order = sortperm(tmp_dist[:,i])
        for j in 1:2^n_bits
            dist[j,i] = tmp_dist[order[j],i]
            matrices[j,i] = all_combs[order[j]]
        end
    end
    write("Matrix_data/fixed_points_eps_"*string(n_bits)*"bits_"
                            *string(epsilon)*"eps_fixedpts", fixed_pts)
    write("Matrix_data/fixed_points_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_distances", dist)
    write("Matrix_data/fixed_points_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_matrices", matrices)
    return fixed_pts, dist, matrices
end


function load_fixed_points(epsilon, n_bits=16)
    fixed_pts = Array{var_type}(3,3,9)
    dist = Array{distance_type}(2^n_bits,9)
    matrices = Array{storage_type}(2^n_bits,9)
    read!("Matrix_data/fixed_points_eps_"*string(n_bits)*"bits_"
                            *string(epsilon)*"eps_fixedpts", fixed_pts)
    read!("Matrix_data/fixed_points_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_distances", dist)
    read!("Matrix_data/fixed_points_eps_"*string(n_bits)*"bits_"
                                *string(epsilon)*"eps_matrices", matrices)
    return fixed_pts, dist, matrices
end


function distance(M1, M2)
    return sqrt(1./2.*sum(angle(eigvals(M1'*M2)).^2))
end


function distances(M1, matrices, gen_mat)
    return map(x->distance(M1,x),
                    [prod([gen_mat[:,:,j]^(div(matrices[i],2^(j-1))%2)
                    for j in 1:size(gen_mat,3)]) for i in 1:length(matrices)])
end


function fix_unitarity!(matrix)
    n_rows, n_columns = size(matrix)
    for n in 1:n_rows, m in 1:n_columns
        if abs(real(matrix[n,m])) < 1e-10
            matrix[n,m] = imag(matrix[n,m])*im
        end
        if abs(imag(matrix[n,m])) < 1e-10
            matrix[n,m] = real(matrix[n,m])+0.0im
        end
    end
    for n in 1:n_rows
        for m in 1:n-1
            matrix[n,:] -= dot(matrix[n,:],matrix[m,:])/
                                     norm(matrix[m,:])*matrix[m,:]
        end
        matrix[n,:] /= norm(matrix[n,:])
    end
    return
end
function fix_unitarity(matrix)
    tmp_matrix = copy(matrix)
    fix_unitarity!(tmp_matrix)
    return tmp_matrix
end


function find_gp_matrices(filename)

    matrices = Array{Array{var_type,2},1}(0)

    for i in gen, j in gen, k in gen, l in gen, m in gen, n in gen, o in gen,
            p in gen, q in gen, r in gen, s in gen
        is_new = true
        mat = i*j*k*l*m*n*o*p*q*r*s
        for a in matrices
            if mat ≈ a
                is_new = false
                break
            end
        end
        if is_new
            push!(matrices, mat)
        end
    end

    storage_array = Array{var_type}(3,3,length(matrices))
    for n in 1:length(matrices)
        storage_array[:,:,n] = matrices[n]
    end

    s = write(filename, storage_array)

    print("Found $(length(matrices)) matrices in the group.")
    println(" Saved in $(filename) ($(s) bytes)")
end


function load_gp_matrices(filename, n_matrices)

    storage_array = Array{var_type}(3,3,n_matrices)
    matrices = Array{Array{var_type,2},1}(n_matrices)

    read!(filename, storage_array)
    for n in 1:n_matrices
        matrices[n] = storage_array[:,:,n]
    end

    return matrices
end


function subdiv_gp_matrices(filename, matrices)

    new_matrices = Array{Array{var_type,2},1}(0)

    n_matrices = length(matrices)
    min_dist = 1.1*min([distance(matrices[1],matrices[i]) for i in 2:n_matrices]...)

    for i in 1:n_matrices, j in i+1:n_matrices
        if distance(matrices[i],matrices[j]) > min_dist
            continue
        end
        # if true#!(distance(matrices[i],matrices[i]*sqrtm(matrices[i]'*matrices[j])) ≈ distance(matrices[i]*sqrtm(matrices[i]'*matrices[j]),matrices[j]))
        #     println(distance(matrices[i],matrices[i]*sqrtm(matrices[i]'*matrices[j])),"  ", distance(matrices[i]*sqrtm(matrices[i]'*matrices[j]),matrices[j]))
        # end
        push!(new_matrices,matrices[i]*sqrtm(matrices[i]'*matrices[j]))
        # tmp_m = logm(matrices[i]'*matrices[j])/im
        # tmp_comps = Array{Float64}(8)
        # tmp_comps[1] = real(tmp_m[2,1])
        # tmp_comps[2] = imag(tmp_m[2,1])
        # tmp_comps[4] = real(tmp_m[3,1])
        # tmp_comps[5] = imag(tmp_m[3,1])
        # tmp_comps[6] = real(tmp_m[3,2])
        # tmp_comps[7] = imag(tmp_m[3,2])
        # tmp_comps[8] = -sqrt(3.)/2.*real(tmp_m[3,3])
        # tmp_comps[3] = real(tmp_m[1,1]) - tmp_comps[8]/sqrt(3.)
        # push!(new_matrices,matrices[i]*expm(im*sum(T_n.*tmp_comps)/2.))
    end

    storage_array = Array{var_type}(3,3,length(new_matrices))
    for n in 1:length(new_matrices)
        storage_array[:,:,n] = new_matrices[n]
    end

    s = write(filename, storage_array)

    print("Found $(length(new_matrices)) matrices in the group.")
    println(" Saved in $(filename) ($(s) bytes)")
end


function multilateration(matrices, fixed_to_all, distance_to_fixed, epsilon=0.2)
    intersection = Array{storage_type}(0)
    n_fixed = length(distance_to_fixed)
    ranges = Array{index_type}(2,n_fixed)
    for i in 1:n_fixed
        ranges[1,i] = searchsortedfirst(fixed_to_all[:,i],
                                                distance_to_fixed[i]-epsilon)
        ranges[2,i] = searchsortedlast(fixed_to_all[:,i],
                                                distance_to_fixed[i]+epsilon)
        if ranges[1,i] == 0 || ranges[1,i] == n_fixed+1 || ranges[2,i] == 0 ||
            ranges[2,i] == n_fixed+1
            #println("One of the intervals is empty")
            return intersection
        end
        #println(ranges[1,i]," ",ranges[2,i],"   ",ranges[2,i]-ranges[1,i])
    end
    range_lengths = [ranges[2,i]-ranges[1,i] for i in 1:n_fixed]
    order = sortperm(range_lengths)

    for i in ranges[1,order[1]]:ranges[2,order[1]]

        mat = matrices[i,order[1]]
        is_in_intersection = true

        for n in order[2:end]
            is_in_range = false
            for j in ranges[1,n]:ranges[2,n]
                if mat == matrices[j,n]
                    is_in_range = true
                    break
                end
            end
            if !is_in_range
                is_in_intersection = false
                break
            end
        end

        if is_in_intersection
            push!(intersection, mat)
        end
    end
    # if length(intersection) == 0
    #     println("Intersection is empty")
    # else
    #     println("Intersection has ",length(intersection)," elements")
    # end
    return intersection#, ranges
end

function multilateration_simd(matrices, fixed_to_all, distance_to_fixed, epsilon = 0.2)
    intersection = Array{storage_type}(0)
    n_fixed = length(distance_to_fixed)
    ranges = Array{index_type}(2,n_fixed)
    for i in 1:n_fixed
        ranges[1,i] = searchsortedfirst(fixed_to_all[:,i],
                                                distance_to_fixed[i]-epsilon)
        ranges[2,i] = searchsortedlast(fixed_to_all[:,i],
                                                distance_to_fixed[i]+epsilon)
        if ranges[1,i] == 0 || ranges[1,i] == n_fixed+1 || ranges[2,i] == 0 ||
            ranges[2,i] == n_fixed+1
            println("One of the intervals is empty")
            return Array{storage_type}(0)
        end
        # println(ranges[2,i]-ranges[1,i])
    end
    range_lengths = [ranges[2,i]-ranges[1,i] for i in 1:n_fixed]
    m, n_small = findmin(range_lengths)

    for i in ranges[1,n_small]:ranges[2,n_small]

        mat = matrices[i,n_small]
        is_in_intersection = true

        for n in 1:n_fixed
            if n == n_small
                continue
            end
            is_in_range = zero(Int32)
            @simd for j in ranges[1,n]:ranges[2,n]
                is_in_range += (mat == matrices[j,n])
            end
            if is_in_range == 0
                is_in_intersection = false
                break
            end
        end

        if is_in_intersection
            push!(intersection, mat)
        end
    end
    # if length(intersection) == 0
    #     println("Intersection is empty")
    #     epsilon += 0.05
    # else
    #     println("Intersection has ",length(intersection)," elements")
    # end
    return intersection
end


function find_closest_matrix_multi(M, gp_matrices, fixed_matrices, all_mat,
                                    dist_to_fixed, gen_mat, eps=0.2)
    gp_i = findmin([distance(M,gp_matrices[j]) for j in 1:1080])[2]
    distf = [distance(M,gp_matrices[gp_i]*fixed_matrices[:,:,j]) for j in 1:9]
    intersection = Array{storage_type}(0)
    ranges = Array{UInt32}(2,9)
    for t in 1:10
        intersection = multilateration(all_mat, dist_to_fixed, distf,eps+0.02t)
        if length(intersection) != 0
            break
        end
    end
    if length(intersection) == 0
        return gp_i, storage_type(0), distance(M,gp_matrices[gp_i]), ranges
    end
    dist, k = findmin(distances(gp_matrices[gp_i]'*M, intersection, gen_mat))
    return gp_i, intersection[k], dist#, ranges
end


function find_closest_matrix_direct(M, gp_matrices, gen_mat)
    gp_i = findmin([distance(M,gp_matrices[j]) for j in 1:1080])[2]
    all_comb = all_bit_combinations(size(gen_mat,3))
    dist, k = findmin(distances(gp_matrices[gp_i]'*M, all_comb, gen_mat))
    return gp_i, all_comb[k], dist
end


end
