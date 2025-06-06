include("../tools/imports.jl")

#Fig. 5: Comparison of OES of related ensembles

f = open("data/perm_MPOs/l15s12.dat", read = true)
d = [parse.(Float64, split(s[2:end-1],",")) for s in split(read(f, String), "\n")[1:end-1]]
close(f)
aut = [vcat(dat, zeros(Float64, 2^12 - length(dat))) for dat in d]
circ1 = aut[1]
aut = vcat(aut...)

f = open("data/automata_vs_bern/automata_spectrum.txt", read = true)
aut = parse.(Float64, split(read(f, String), "\n")[2:end-1])
aut1 = aut[1:2^12]
close(f)

f = open("data/automata_vs_bern/bernoulli_spectrum.txt", read = true)
bern = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

f = open("data/automata_vs_bern/permutation_spectrum.txt", read = true)
perm = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

bwidth = .03175
blist = 0:bwidth:3
histogram(bern, bins = blist, normalize = true, label = "Bernoulli avg OES", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(aut, bins = blist, normalize = true, label = L"$\{P^tXP\}_P$ avg OES", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(perm, bins = blist, normalize = true, label = L"$\{P\}$ avg OES", linetype = :stephist, linewidth = 2)
histogram!(aut1, bins = blist, normalize = true, label = L"$P^tXP$ OES", linetype = :stephist, linewidth = 2)
est = exp(-exp(-1)) + 2*exp(-1) -1
hline!([est/bwidth], label = L"$p(\{0\})$ pred.", color = :black, linestyle= :dash, linewidth = 3)
xlims!(0,2.2)
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
savefig("final_paper_figures/aut_vs_bern_vs_perm.svg")

#Inset of Figs. 7 and 8
bwidth = .03175
blist = 0:bwidth:3
p = plot(size = (500,200), axesfontsize = 12, guidefontsize=12)
histogram!(aut1, bins = blist, normalize = true, label = L"N=12, \text{one sample}", linetype = :stephist, palette = seaborn, linewidth = 2)
#histogram!(aut, bins = blist, normalize = true, label = L"N=12, \text{one sample}", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(bern, bins = blist, normalize = true, label = "Bernoulli spectrum", linetype = :stephist, palette = seaborn, linewidth = 2)
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
yticks!(0:2:12)
savefig("final_paper_figures/aut_spectrum_inset.svg")
savefig("final_paper_figures/aut_circuit_vs_bern_insert.svg")