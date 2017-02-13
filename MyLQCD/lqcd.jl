module lqcd
#using ProgressMeter

# Set SU(n) group and varable type
const su_n = 3
const var_type = Complex{Float64}
const id = eye(var_type,su_n)
const ze = zeros(id)
const epsilon = 0.18

# Pauli matrices
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

export new_gauge_field, update, measure, save_gauge_field, load_gauge_field,
       save_comp_gauge_field, load_comp_gauge_field

function new_gauge_field(Ns, Nt, mode="cold")

    field = Array{Array{var_type,2}}(4,Ns,Ns,Ns,Nt)

    if mode == "cold"
        for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
            field[d,x,y,z,t] = copy(id)
        end
    elseif mode == "hot"
        if su_n == 2
            for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
                r = rand(3)
                field[d,x,y,z,t] = expm(im*(r[1]*sigma_1+r[2]*sigma_2
                                            +r[3]*sigma_3))
            end
        elseif su_n == 3
            for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
                r = rand(8)
                field[d,x,y,z,t] = expm((r[1]*T_1+r[2]*T_2+r[3]*T_3+r[4]*T_4
                                         +r[5]*T_5+r[6]*T_6+r[7]*T_7+r[8]*T_8)
                                         *im)
            end
        else
            println("Hot mode not implemented for other groups")
        end
    end

    return field
end


function random_matrix(A, beta)

    if su_n == 3
        if size(A,1) != 2 || size(A,2) != 2
            println("Input must be 2x2 block")
        end

        a = sqrt(abs(det(A)))
        V = A/a

        x0 = 0.
        x = [0.,0.,0.]

        while true
            r = rand()
            r1 = 1. - rand()
            r2 = 1. - rand()
            r3 = 1. - rand()

            lambda2 = -1./(2.*a*beta)*(log(r1) + cos(2.0pi*r2)^2*log(r3))

            if r^2 > 1. - lambda2
                continue
            end

            x0 = 1. - 2.0lambda2
            break
        end

        while true
            x = [2.0*rand()-1,2.0*rand()-1,2.0*rand()-1]

            if norm(x) > 1
                continue
            end

            x *= sqrt(1.-x0^2)/norm(x)
            break
        end

        X = x0*sigma_0 + x[1]*sigma_1*im + x[2]*sigma_2*im + x[3]*sigma_3*im

    elseif su_n == 2
        println("SU(2) not yet implemented")
    else
        println("Other gauge groups not implemented")
    end

    return X
end
function random_matrix()
    r0 = rand()-1.
    r = [rand()-1.,rand()-1.,rand()-1.]
    x = epsilon*r/norm(r)
    x0 = sign(r0)sqrt(1.-epsilon^2)
    X = x0*sigma_0 + x[1]*sigma_1*im + x[2]*sigma_2*im + x[3]*sigma_3*im
    return X
end


function single_update(field, x, y, z, t, d, Ns, Nt, beta, mode)

    # Sum of staples
    A = copy(ze)
    updated = false

    for dir in 1:4
        if d == dir
            continue
        end

        tmp_x, tmp_y, tmp_z, tmp_t = x, y, z, t

        if d == 1
            tmp_x = mod(tmp_x,Ns)+1
        elseif d == 2
            tmp_y = mod(tmp_y,Ns)+1
        elseif d == 3
            tmp_z = mod(tmp_z,Ns)+1
        elseif d == 4
            tmp_t = mod(tmp_t,Nt)+1
        end

        tmp_A1 = field[dir, tmp_x, tmp_y, tmp_z, tmp_t]

        if dir == 1
            tmp_x = mod(tmp_x,Ns)+1
        elseif dir == 2
            tmp_y = mod(tmp_y,Ns)+1
        elseif dir == 3
            tmp_z = mod(tmp_z,Ns)+1
        elseif dir == 4
            tmp_t = mod(tmp_t,Nt)+1
        end
        if d == 1
            tmp_x = mod(tmp_x-2,Ns)+1
        elseif d == 2
            tmp_y = mod(tmp_y-2,Ns)+1
        elseif d == 3
            tmp_z = mod(tmp_z-2,Ns)+1
        elseif d == 4
            tmp_t = mod(tmp_t-2,Nt)+1
        end

        tmp_A1 *= field[d, tmp_x, tmp_y, tmp_z, tmp_t]'

        if dir == 1
            tmp_x = mod(tmp_x-2,Ns)+1
        elseif dir == 2
            tmp_y = mod(tmp_y-2,Ns)+1
        elseif dir == 3
            tmp_z = mod(tmp_z-2,Ns)+1
        elseif dir == 4
            tmp_t = mod(tmp_t-2,Nt)+1
        end

        tmp_A1 *= field[dir, tmp_x, tmp_y, tmp_z, tmp_t]'

        if d == 1
            tmp_x = mod(tmp_x,Ns)+1
        elseif d == 2
            tmp_y = mod(tmp_y,Ns)+1
        elseif d == 3
            tmp_z = mod(tmp_z,Ns)+1
        elseif d == 4
            tmp_t = mod(tmp_t,Nt)+1
        end
        if dir == 1
            tmp_x = mod(tmp_x-2,Ns)+1
        elseif dir == 2
            tmp_y = mod(tmp_y-2,Ns)+1
        elseif dir == 3
            tmp_z = mod(tmp_z-2,Ns)+1
        elseif dir == 4
            tmp_t = mod(tmp_t-2,Nt)+1
        end

        tmp_A2 = field[dir, tmp_x, tmp_y, tmp_z, tmp_t]'

        if d == 1
            tmp_x = mod(tmp_x-2,Ns)+1
        elseif d == 2
            tmp_y = mod(tmp_y-2,Ns)+1
        elseif d == 3
            tmp_z = mod(tmp_z-2,Ns)+1
        elseif d == 4
            tmp_t = mod(tmp_t-2,Nt)+1
        end

        tmp_A2 *= field[d, tmp_x, tmp_y, tmp_z, tmp_t]'
        tmp_A2 *= field[dir, tmp_x, tmp_y, tmp_z, tmp_t]

        A += tmp_A1 + tmp_A2
    end

    if mode == "heatbath"
        W = field[d,x,y,z,t] * A
        R = copy(id)
        R[[1,2],[1,2]] = random_matrix(W[[1,2],[1,2]], beta)

        W = R*W
        S = copy(id)
        S[[1,3],[1,3]] = random_matrix(W[[1,3],[1,3]], beta)

        W = S*W
        T = copy(id)
        T[[2,3],[2,3]] = random_matrix(W[[2,3],[2,3]], beta)

        field[d,x,y,z,t] = T*S*R*field[d,x,y,z,t]
        updated = true
    elseif mode == "metropolis"
        R = copy(id)
        S = copy(id)
        T = copy(id)
        R[[1,2],[1,2]] = random_matrix()
        S[[1,3],[1,3]] = random_matrix()
        T[[2,3],[2,3]] = random_matrix()

        X = R*S*T
        X = (rand() < 0.5 ? X : X')

        deltaS = -beta/su_n*real(trace((X*field[d,x,y,z,t]-field[d,x,y,z,t])*A))

        r = rand()
        if r <= exp(-deltaS)
            field[d,x,y,z,t] = X*field[d,x,y,z,t]
            updated = true
        end
    end
    return updated
end


function update(field, beta, n_sweeps, mode="metropolis")
    Ns = size(field,2)
    Nt = size(field,5)

    accepted = 0

    for n in 1:n_sweeps, t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
        if single_update(field, x, y, z, t, d, Ns, Nt, beta, mode)
            accepted += 1
        end
    end
    println("Accepted ratio: ",accepted/(n_sweeps*Nt*Ns^3*4))
end


function measure(field)
    Ns = size(field,2)
    Nt = size(field,5)

    av_plaq, av_plaq_s, av_plaq_t = 0., 0., 0.
    n_plaq, n_plaq_s, n_plaq_t = 0, 0, 0

    for dir1 in 1:4, dir2 in (dir1+1):4, t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns

        tmp_x, tmp_y, tmp_z, tmp_t = x, y, z, t

        plaq = field[dir1, tmp_x, tmp_y, tmp_z, tmp_t]

        if dir1 == 1
            tmp_x = mod(tmp_x,Ns)+1
        elseif dir1 == 2
            tmp_y = mod(tmp_y,Ns)+1
        elseif dir1 == 3
            tmp_z = mod(tmp_z,Ns)+1
        elseif dir1 == 4
            tmp_t = mod(tmp_t,Nt)+1
        end

        plaq *= field[dir2, tmp_x, tmp_y, tmp_z, tmp_t]

        if dir2 == 1
            tmp_x = mod(tmp_x,Ns)+1
        elseif dir2 == 2
            tmp_y = mod(tmp_y,Ns)+1
        elseif dir2 == 3
            tmp_z = mod(tmp_z,Ns)+1
        elseif dir2 == 4
            tmp_t = mod(tmp_t,Nt)+1
        end
        if dir1 == 1
            tmp_x = mod(tmp_x-2,Ns)+1
        elseif dir1 == 2
            tmp_y = mod(tmp_y-2,Ns)+1
        elseif dir1 == 3
            tmp_z = mod(tmp_z-2,Ns)+1
        elseif dir1 == 4
            tmp_t = mod(tmp_t-2,Nt)+1
        end

        plaq *= field[dir1, tmp_x, tmp_y, tmp_z, tmp_t]'

        if dir2 == 1
            tmp_x = mod(tmp_x-2,Ns)+1
        elseif dir2 == 2
            tmp_y = mod(tmp_y-2,Ns)+1
        elseif dir2 == 3
            tmp_z = mod(tmp_z-2,Ns)+1
        elseif dir2 == 4
            tmp_t = mod(tmp_t-2,Nt)+1
        end

        plaq *= field[dir2, tmp_x, tmp_y, tmp_z, tmp_t]'

        tr = 1./3.*real(trace(plaq))

        av_plaq += tr
        n_plaq += 1
        if dir1 != 4 && dir2 != 4
            av_plaq_s += tr
            n_plaq_s += 1
        elseif dir1 == 4 || dir2 == 4
            av_plaq_t += tr
            n_plaq_t += 1
        end
    end

    av_plaq /= n_plaq
    av_plaq_s /= n_plaq_s
    av_plaq_t /= n_plaq_t

    println("Average: ",av_plaq)
    println("Average spatial: ",av_plaq_s)
    println("Average temporal: ",av_plaq_t,"\n")

    return av_plaq
end


function save_gauge_field(field, filename)
    Ns = size(field,2)
    Nt = size(field,5)
    tmp_field = Array{var_type}(su_n,su_n,4,Ns,Ns,Ns,Nt)
    for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
        tmp_field[:,:,d,x,y,z,t] = field[d,x,y,z,t]
    end
    s = write(filename, tmp_field)
    println("Saved gauge field in $(filename) ($(s) bytes)")
end


function load_gauge_field(filename, Ns, Nt)
    field = Array{Array{var_type,2}}(4,Ns,Ns,Ns,Nt)
    tmp_field = Array{var_type}(su_n,su_n,4,Ns,Ns,Ns,Nt)
    read!(filename, tmp_field)
    for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
        field[d,x,y,z,t] = tmp_field[:,:,d,x,y,z,t]
    end
    println("Loaded gauge field from $(filename)")
    return field
end


function save_comp_gauge_field(field, filename)
    Ns = size(field,2)
    Nt = size(field,5)
    tmp_comp_field = Array{real(var_type)}(8,4,Ns,Ns,Ns,Nt)
    for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
        tmp_m = logm(field[d,x,y,z,t])/im
        tmp_comp_field[1,d,x,y,z,t] = real(tmp_m[2,1])
        tmp_comp_field[2,d,x,y,z,t] = imag(tmp_m[2,1])
        tmp_comp_field[4,d,x,y,z,t] = real(tmp_m[3,1])
        tmp_comp_field[5,d,x,y,z,t] = imag(tmp_m[3,1])
        tmp_comp_field[6,d,x,y,z,t] = real(tmp_m[3,2])
        tmp_comp_field[7,d,x,y,z,t] = imag(tmp_m[3,2])
        tmp_comp_field[8,d,x,y,z,t] = -sqrt(3.)/2.*real(tmp_m[3,3])
        tmp_comp_field[3,d,x,y,z,t] = real(tmp_m[1,1]) - (
                                           tmp_comp_field[8,d,x,y,z,t])/sqrt(3.)
    end
    s = write(filename, tmp_comp_field)
    println("Saved gauge field in $(filename) ($(s) bytes)")
end


function load_comp_gauge_field(filename, Ns, Nt)
    field = Array{Array{var_type,2}}(4,Ns,Ns,Ns,Nt)
    tmp_comp_field = Array{real(var_type)}(8,4,Ns,Ns,Ns,Nt)
    read!(filename, tmp_comp_field)
    for t in 1:Nt, z in 1:Ns, y in 1:Ns, x in 1:Ns, d in 1:4
        field[d,x,y,z,t] = expm(im*sum(tmp_comp_field[:,d,x,y,z,t].*T_n))
    end
    println("Loaded gauge field from $(filename)")
    return field
end


end
