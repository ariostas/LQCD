push!(LOAD_PATH, pwd())
using lqcd

Ns = 4
Nt = 10
beta = 6.0
n_therm = 10
n_corr = 1
n_meas = 100

field = new_gauge_field(Ns, Nt)

update(field, beta, n_therm)

for n in 1:n_meas
    measure(field)
    if n < n_meas
        update(field, beta, n_corr)
    end
end

save_gauge_field(field, "testsave.qcd")

field2 = load_gauge_field("testsave.qcd",Ns,Nt)

measure(field2)

save_comp_gauge_field(field, "testsave2.qcd")

field3 = load_comp_gauge_field("testsave2.qcd",Ns,Nt)

measure(field3)
