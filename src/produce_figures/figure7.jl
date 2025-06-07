include("../tools/imports.jl")

#Plot Fig. 7
f = open("data/rand_perm_moms.dat", "r")
moms = Dict()
sizes = [6,8,10,12]
ks = [2,3,4,5]
for (size,s) in zip(sizes,split(read(f, String), "\n")[1:end-1])
    arr = [t[1:end] for t in split.(s[14:end-2], "[")]
    moms[size] = [parse.(Float64,strip.(s, ']')[1:end-1]) for s in split.(arr, ",")]
end

#Data for a permutation on system size N=14 (computed separately due to memory constraints)
f = open("data/X_data_14.dat", "r")

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

p = plot(guidefontsize = 14)
moms_ideal = [3, 12, 57, 303]
for (k,m) in zip(1:4, moms_ideal)
    plot!(sizes, [(m-mean(moms[s][k]))/m for s in sizes], yerr = [std(moms[s][k])/(m*sqrt(length(moms[s][k]))) for s in sizes], yaxis = :log, label = "k = $(k+1)", color = seaborn[k], linewidth = 2, marker = "x")
end

xlabel!(L"N")
ylabel!(L"\frac{\langle \lambda^k \rangle - m_k}{m_k}")
title!(L"Convergence of $X(\tau)$ OES moments to Bernoulli ensemble")

savefig("figures/fig7.svg")

f = open("data/automata_vs_bern/automata_spectrum.txt", read = true)
aut = parse.(Float64, split(read(f, String), "\n")[2:end-1])
aut1 = aut[1:2^12]
close(f)

f = open("data/automata_vs_bern/bernoulli_spectrum.txt", read = true)
bern = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

#Inset of Fig. 7
bwidth = .03175
blist = 0:bwidth:3
p = plot(size = (500,200), axesfontsize = 12, guidefontsize=12)
histogram!(aut1, bins = blist, normalize = true, label = L"N=12, \text{one sample}", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(bern, bins = blist, normalize = true, label = "Bernoulli spectrum", linetype = :stephist, palette = seaborn, linewidth = 2)
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
yticks!(0:2:12)
savefig("figures/fig7_inset.svg")
