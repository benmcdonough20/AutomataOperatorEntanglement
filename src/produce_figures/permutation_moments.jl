include("../tools/imports.jl")

#Plot Fig. 6
f = open("data/off_diagonal/diagonal.dat")
dat = read(f, String)
close(f)
dat1 = split(dat, "\n")[1:end-1]
dist = parse.(Float64, dat1)

histogram(dist./8, normalize = true, linetype = :stephist, label = L"2^{-N/4} \times P^tZP", palette = seaborn, bins = 0:.1:2.1)
plot!(LinRange(0, 2, 100), x->1/π * sqrt(4-x^2)/64, label = "semicircle (rescaled)", color = palette(:seaborn_deep)[4])
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
ylims!(0, .04)
yticks!(0:.004:02)
savefig("figures/Z_aut_evo.svg")

#(Not included in paper) Show that OES of automaton-evolved OES is the same for all off-diagonal operators
f = open("data/off_diagonal/spectrum_vs_angle.dat", read = true)
dat = read(f, String)
close(f)
dists = [parse.(Float64, split(d[2:end-1], ",")) for d in split(dat, "\n")[1:end-1]]

c1 = HSV(palette(:seaborn_deep)[1])
cols = [HSV(c1.h, c1.v, c1.s -.15*i) for i in -2:2]
p = plot(size = (300,400), palette = :seaborn_deep)
for (i,d) in enumerate(dists)
    histogram!(d,normalize = true, linetype = :stephist, color = cols[i])
end
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
savefig("final_paper_figures/off_diagonal.svg")

#Plot Fig. 7
f = open("rand_perm_moms.dat", "w")
println(f, moms)
close(f)


for s in sizes
    for i in 1:4
        moms[s][i] = real(moms[s][i])
    end
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
display(p)

xlabel!(L"N")
ylabel!(L"\frac{\langle \lambda^k \rangle - m_k}{m_k}")
title!(L"Convergence of $X(\tau)$ OES moments to Bernoulli ensemble")

savefig("figs/paper/moments_convergence.svg")