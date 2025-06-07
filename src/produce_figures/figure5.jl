include("../tools/imports.jl")

#Fig. 5: Comparison of OES of related ensembles

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
savefig("figures/fig5.svg")


