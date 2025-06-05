using Pkg; Pkg.activate("automata")
using LaTeXStrings

include("../tools/plotting_defaults.jl")

using Plots

f = open("data/perm_MPOs/l15s12.dat", read = true)
d = [parse.(Float64, split(s[2:end-1],",")) for s in split(read(f, String), "\n")[1:end-1]]
close(f)
aut = [vcat(dat, zeros(Float64, 2^12 - length(dat))) for dat in d]
aut1 = aut[1]
aut = vcat(aut...)

f = open("data/automata_vs_bern/automata_spectrum.txt", read = true)
aut = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

f = open("data/automata_vs_bern/bernoulli_spectrum.txt", read = true)
bern = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

f = open("data/automata_vs_bern/permutation_spectrum.txt", read = true)
perm = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

bwidth = .03175
bins = 0:bwidth:3
histogram(bern, bins = bins, normalize = true, label = "Bernoulli avg OES", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(aut, bins = bins, normalize = true, label = L"$\{P^tXP\}_P$ avg OES", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(perm, bins = bins, normalize = true, label = L"$\{P\}$ avg OES", linetype = :stephist, linewidth = 2)
histogram!(aut1, bins = bins, normalize = true, label = L"$P^tXP$ OES", linetype = :stephist, linewidth = 2)
est = exp(-exp(-1)) + 2*exp(-1) -1
hline!([est/bwidth], label = L"$p(\{0\})$ pred.")
plot!(LinRange(0,2,100), x->sqrt(1/pi)*sqrt(4-x^2))
xlims!(0,2.2)
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
savefig("final_paper_figures/aut_vs_bern_vs_perm.svg")


plot(size=(500,200), palette = :seaborn_dark)
histogram!(aut1, bins = bins, normalize = true, label = "N=12, one sample", linetype = :stephist, linewidth = 2)
pred = (exp(-exp(-1)) + 2*exp(-1)-1)/bwidth
hline!([pred], label = "predicted kernel")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
yticks!(0:2:12)
savefig("final_paper_figures/aut_spectrum_inset.svg")

println(sum(aut1 .< 1E-13)/length(aut1))