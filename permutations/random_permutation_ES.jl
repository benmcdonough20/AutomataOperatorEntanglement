using Pkg; Pkg.activate("automata")
using ProgressBars

include("../tools/plotting_defaults.jl")
include("../tools/permutation_tools.jl")
include("../tools/ES_tools.jl")

using LaTeXStrings
using LinearAlgebra

#average permutation spectrum Bernoulli vs Automaton
samples = 2^6
ss = 12
Px = [0.0 1.0; 1.0 0.0]
d = 2^ss
X = kron([[[1 0; 0 1] for i in 1:Int(ss/2-1)];[Px];[[1 0; 0 1] for i in 1:Int(ss/2)]]...)

dist_bern = []
for i in ProgressBar(1:samples)
    P = (rand(d,d) .<= 1/d).*exp.(rand([1im,-1im], (d,d)))
    ptranspose!(P)
    push!(dist_bern, svdvals(P)...)
end

dist_aut = []
for i in ProgressBar(1:samples)
    P = rand_perm_mat(d);
    O = P*X*P'
    ptranspose!(O)
    push!(dist_aut, svdvals(O)...)
end

dist_perm = []
for i in ProgressBar(1:samples)
    P = rand_perm_mat(d);
    ptranspose!(P)
    push!(dist_perm, svdvals(P)...)
end

#save data
f = open("data/automata_vs_bern/bernoulli_spectrum.txt", write = true, create = true)
println(f, "samples: 64, system size: 12, matrix type: bernoulli")
for d in dist_bern
println(f, d)
end
close(f)

f = open("data/automata_vs_bern/automata_spectrum.txt", write = true, create = true)
println(f, "samples: 64, system size: 12, matrix type: automata")
for d in dist_aut
println(f, d)
end
close(f)

f = open("data/automata_vs_bern/permutation_spectrum.txt", write = true, create = true)
println(f, "samples: 64, system size: 12, matrix type: permutation")
for d in dist_perm
println(f, d)
end
close(f)

D = 12
ops = [kron([0 exp(1im * π/8*i); exp(-1im * π/8*i) 0],Matrix(I, 2^(D-1), 2^(D-1))) for i in 0:4]

dists = [[] for op in ops]
samples = 64
for (j,op) in enumerate(ops)
    for i in ProgressBar(1:samples)
        P = rand_perm_mat(2^D);
        M = P*op*P'
        ptranspose!(M)
        push!(dists[j], svdvals(M)...)
    end
end
dists
f = open("data/off_diagonal/spectrum_vs_angle.dat", write = true)
for d in dists
    println(f, convert.(Float64,d))
end
close(f)

f = open("data/off_diagonal/spectrum_vs_angle.dat", read = true)
dat = read(f, String)
close(f)
dists = [parse.(Float64, split(d[2:end-1], ",")) for d in split(dat, "\n")[1:end-1]]

c1 = HSV(palette(:seaborn_deep)[1])
cols = [HSV(c1.h, c1.v, c1.s -.15*i) for i in -2:2]
p = plot(size = (300,400), palette = :seaborn_deep)
for (i,d) in enumerate(dists)
    histogram!(d, bins = bins, normalize = true, linetype = :stephist, color = cols[i])
end
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
savefig("final_paper_figures/off_diagonal.svg")

op = kron([1 0; 0 -1],Matrix(I, 2^(D-1), 2^(D-1)))
dist = []
samples = 64
for i in ProgressBar(1:samples)
    P = rand_perm_mat(2^D);
    M = P*op*P'
    ptranspose!(M)
    push!(dist, svdvals(M)...)
end

f = open("data/off_diagonal/diagonal.dat", write = true)
for d in dist
    println(f, d)
end
close(f)

f = open("data/off_diagonal/diagonal.dat")
dat = read(f, String)
close(f)
dat1 = split(dat, "\n")[1:end-1]
dist = parse.(Float64, dat1)

#histogram(dist[dist .> 1E-13]./8, normalize = true, linetype = :stephist, label = L"2^{-N/4} \times P^tZP", palette = seaborn, bins = 0:.05:2.1)
histogram(dist./8, normalize = true, linetype = :stephist, label = L"2^{-N/4} \times P^tZP", palette = seaborn, bins = 0:.1:2.1)
plot!(LinRange(0, 2, 100), x->1/π * sqrt(4-x^2)/64, label = "semicircle (rescaled)", color = palette(:seaborn_deep)[4])
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
ylims!(0, .04)
yticks!(0:.004:02)

savefig("final_paper_figures/Z_aut_evo_unscaled.svg")
savefig("final_paper_figures/Z_aut_evo_scaled.svg")

savefig("final_paper_figures/Z_aut_evo.svg")

sum(dist .> 1E-5)/length(dist)

D = 10
op = kron([1 0; 0 -1],Matrix(I, 2^(D-1), 2^(D-1)))
H = kron([1/sqrt(2)*[1 1; 1 -1] for i in 1:D]...)
P = rand_perm_mat(2^D)
M = P*op*P'
ptranspose!(M)
ns = [norm(M[i,:])^2 for i in 1:2^D]
plot(ns)
M

lst = svdvals(M)
histogram(lst[lst .> 1E-10])

M = randn((2^10, 2^10)) .* 1/2^5
M = diagm(randn(2^10))
ptranspose!(M)
v = M[1,:]
norm(v)^2
lst = svdvals(M)
D = 10
ns = [norm(M[i,:])^2 for i in 1:2^D]
plot(ns)
lst




samples = 2^16
ss = [6,8,10,12]
atoms = []
stds = []
for s in ss
    N = 2^s
    dist_perm = []
    for i in ProgressBar(1:samples * 2.0^(-s))
        #P = rand_perm_mat(N);
        #Z,M = gen_large_cluster(2^s)
        M = rand(2^s, 2^s) .<= 1/2^s
        push!(dist_perm, svdvals(M))
    end
    ests = []
    for i in 1:length(dist_perm)
        samp = dist_perm[i]
        mp = countmap(round.(samp, digits = 10))
        push!(ests, sum([mp[k] for k in collect(keys(mp)) if mp[k]>2])/length(samp))
    end
    mp = countmap(round.(vcat(dist_perm...), digits = 10))
    push!(atoms, sum([mp[k] for k in collect(keys(mp)) if mp[k]>2])/length(vcat(dist_perm...)))
    push!(stds, std(ests)/sqrt(length(dist_perm)))
end

plot(6:2:12, atoms, linewidth = 2, yerr = stds)
title!("Size of singular component of spectrum")
xlabel!("N")
ylabel!("p")

plot!(LinRange(6, 12, 100), x->1/sqrt((x-6)))
savefig("figs/paper/singularfrac.svg")

#Compute and compare moments
moms = Dict()
sizes = [6,8,10,12, 14]
for s in sizes
    moms[s] = [[],[],[],[]]
end

samples = 2^20 .* 2.0 .^ (-sizes)
for (k,s) in ProgressBar(enumerate(sizes))
    X = kronecker(Matrix(I, 2^(s÷2-1), 2^(s÷2-1)), [0 1; 1 0], Matrix(I, 2^(s÷2), 2^(s÷2)))
    for i in ProgressBar(1:samples[k])
        P = rand_perm_mat(2^s)
        P = P'*X*P
        ptranspose!(P)
        B = P*P'
        M = B
        for j in 1:4
            M = M*B
            push!(moms[s][j],1/2^s*tr(M))
        end
    end
end

f = open("rand_perm_moms.dat", "w")
println(f, moms)
close(f)

for s in sizes
    for i in 1:4
        moms[s][i] = real(moms[s][i])
    end
end

f = open("MPOs/X_data_14.dat", "r")

d = read(f, String)
arrs_str = split(d, "\n")
arrs_str = [a[2:end-1] for a in arrs_str]
arrs = split.(arrs_str, ",")
arrs = [parse.(Float64, a) for a in arrs[1:end-1]]
close(f)

moms[14] = [[],[],[], []]
for a in arrs
    for i in 1:4
        push!(moms[14][i], sum(a .^ ((i+1)*2)) / 2^14)
    end
end

include("tools/plotting_defaults.jl")

p = plot(guidefontsize = 14)
moms_ideal = [3, 12, 57, 303]
for (k,m) in zip(1:4, moms_ideal)
    plot!(sizes, [(m-mean(moms[s][k]))/m for s in sizes], yerr = [std(moms[s][k])/(m*sqrt(length(moms[s][k]))) for s in sizes], yaxis = :log, label = "k = $(k+1)", color = seaborn[k], linewidth = 2, marker = "x")
end
display(p)

xlabel!(L"N")
ylabel!(L"\frac{\langle \lambda^k \rangle - m_k}{m_k}")
title!(L"Convergence of $X(\tau)$ OES moments to Bernoulli ensemble")

savefig("figs/paper/moments_convergence.svg")

log10((303-mean(moms[14][4]))/303)