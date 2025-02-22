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

f = open("data/automata_vs_bern/bernoulli_spectrum.txt", read = true)
bern = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

bwidth = .0375
bins = 0:bwidth:3
plot(size=(500,200), palette = :seaborn_dark)
histogram!(aut, bins = bins, normalize = true, label = "N=12, l=15", linetype = :stephist, alpha = .5, linewidth = 2)
histogram!(bern, bins = bins, normalize = true, label = "Bernoulli spectrum", linetype = :stephist, alpha = .5, linewidth = 2)
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
savefig("final_paper_figures/aut_vs_bern_inset.svg")


plot(size=(500,200), palette = :seaborn_dark)
histogram!(aut1, bins = bins, normalize = true, label = "N=12, one sample", linetype = :stephist, linewidth = 2)
pred = (exp(-exp(-1)) + 2*exp(-1)-1)/bwidth
hline!([pred], label = "predicted kernel")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
yticks!(0:2:12)
savefig("final_paper_figures/aut_spectrum_inset.svg")

println(sum(aut1 .< 1E-13)/length(aut1))