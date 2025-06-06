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