using Pkg; Pkg.activate("automata")
using StatsBase
using Plots
using LaTeXStrings
using LsqFit
using CurveFit
using ProgressBars
using Colors

wigner_PDF = x-> 0 <= x <= 2 ? 1/π * sqrt(4-x^2) : 0 #Integration convention: use right sums!
wigner_CDF = x-> 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)

function read_dist(N, nh)
	f = open("semicircle_data/X_N$(N)_nH$(nh)", read = true)
	l = read(f, String)
	close(f)
	dists_str = split(l,"sample")[2:end]
	dists_num = split.(dists_str, "\n")
	dists = getindex.(dists_num, (2:length(dists_num[1])-1,))
	map(x->parse.(Float64, x), dists)
end

#Wass distance

#note: sample size makes an important difference in this metric, but bin size does not
function get_dist_wass(dist)
	ECDF = ecdf(dist)
	step = 2/(2*length(dist)^(1/3))
	maxstep = ceil(maximum(dist)/step)
	sum([abs(ECDF(i*step) - wigner_CDF(i*step)) for i in 1:maxstep]) * step
end

sc = palette([:red, :blue], 4)

p = plot(yaxis = :log, dpi = 250)

y = [get_dist_wass(vcat(read_dist(N,round(Int, N*2/3))...)) for N in 6:2:12]
b1,b2 = linear_fit(6:2:12, log.(y))

y = [get_dist_wass(vcat(read_dist(12,nh)...)) for nh in 1:4]
c1,c2 = linear_fit(1:4, log.(y))

for n in ProgressBar(6:2:12)
	y = [vcat(read_dist(n, nH)...) for nH in 1:min(n, 11)]
	xax = collect(1:min(n,11))
	x = xax .- (b2*n + b1)/c2
	plot!(x, [get_dist_wass(dist) for dist in y] ./ exp(b2*n+b1), color = sc[(n-4)÷2], label = "N=$(n)")
end

title!("Universal scaling law (X)")
xlabel!(L"n_H - (.667N+1.1)/1.37")
ylabel!(L"W1/exp(-0.667N-1.1)")
display(p)
#aim for 1 million samples!
display(p)

function get_dist_DKL(N, n)
	dist = vcat(read_dist(N,n)...)
	binnum = floor(length(dist)^(1/3)*2)
	binsize = maximum([maximum(dist),2])/binnum
	P = [sum(0 .<= dist .- binsize*i .< binsize) for i in 0:binnum-1] ./ length(dist)
	Q = [wigner_CDF((i+1)*binsize) - wigner_CDF(i*binsize) for i in 0:binnum-1]
	#Q = [wigner_PDF((i+1)*binsize)*binsize for i in 0:binnum] #doesn't work as well for really small values...
	arr = Q .* log.(Q ./ P)
	arr[any.(isnan, arr)] .= 0
	arr[any.(isinf, arr)] .= 0
	return sum(arr)
end

p = plot(yaxis = :log)

for n in ProgressBar(6:2:12)
	plot!([get_dist_DKL(n, nH) for nH in 1:min(n,11)])
end

display(p)

function get_dist_KS(N, n)
	dist = sort(vcat(read_dist(N, n)...))
	ECDF = ecdf(dist)
	maximum([abs(ECDF(dist[x]) - wigner_CDF(dist[x])) for x in 1:length(dist)])
end

p = plot(yaxis = :log)

for n in 6:2:10
	plot!([get_dist_KS(n, nH) for nH in 1:n])
end

display(p)

#https://www.netlib.org/lapack/lug/node96.html precision: 1.3 * 10^-6
#atoms, if they exist, are hard to find--better to work with ranges. For instance, plot the kernel range?

#level spacings
spacings(dist) = (d = sort(dist);[(d[i+2]-d[i+1])/(d[i+1]-d[i]) for i in 1:length(dist)-2])


function linear_cross_entropy(N, n) #pretty nice
	dist = vcat(read_dist(N,n)...)
	return abs(1 - 3π^2/16 *mean([wigner_PDF(x) for x in dist]))
end

function linear_cross_entropy(dist) #pretty nice
	return abs(1 - 3π^2/16 *mean([wigner_PDF(x) for x in dist]))
end

p = plot(yaxis = :log)

for n in 6:2:12
	plot!([linear_cross_entropy(n, nH) for nH in 1:min(n,11)])
end

p = plot(yaxis = :log, dpi = 250)

y = [linear_cross_entropy(vcat(read_dist(N,round(Int, N*2/3))...)) for N in 6:2:12]
b1,b2 = linear_fit(6:2:12, log.(y))

y = [linear_cross_entropy(vcat(read_dist(12,nh)...)) for nh in 1:4]
c1,c2 = linear_fit(1:4, log.(y))

for n in ProgressBar(6:2:12)
	y = [vcat(read_dist(n, nH)...) for nH in 1:min(n, 11)]
	xax = collect(1:min(n,11))
	x = xax .- (b2*n + b1)/c2
	plot!(x, [linear_cross_entropy(dist) for dist in y] ./ exp(b2*n+b1), color = sc[(n-4)÷2], label = "N=$(n)")
end

display(p)

function peak_at_zero(N, n)
	dist = vcat(read_dist(N,n)...)
	dist = round.(dist, digits = 6)
	cutoff = length(dist)^(-1/3)
	return abs(sum(dist .<= cutoff)/length(dist) - wigner_CDF(cutoff))
end

p = plot(yaxis = :log)
for n in 6:2:12
	plot!([peak_at_zero(n, nH) for nH in 1:min(n, 11)])
end

display(p)

plot([peak_at_zero(n, 6) for n in 6:2:10], yaxis = :log)

function tail_prob(N, n) #good
	dist = vcat(read_dist(N,n)...)
	return sum(dist .>= 2)/length(dist)
end

p = plot(yaxis = :log)
for n in 6:2:12
	plot!([tail_prob(n, nH) for nH in 1:min(n, 11)])
end

display(p)

histogram(vcat(read_dist(10,1)...), normalize = true, label = "", dpi = 250)
title!(L"$n_H = 1$")
xlabel!(L"\lambda")
xlims!(0, 2.5)
ylabel!("p")
savefig("figs/more_samples/1.png")

histogram(vcat(read_dist(10,2)...), normalize = true, label = "", dpi = 250)
title!(L"$n_H = 2$")
xlabel!(L"\lambda")
xlims!(0, 2.5)
ylabel!("p")
savefig("figs/more_samples/2.png")


histogram(vcat(read_dist(10,3)...), label = "", dpi = 250, normalize = true)
title!(L"$n_H = 3$")
xlabel!(L"\lambda")
xlims!(0, 2.5)
ylabel!("p")
savefig("figs/more_samples/3.png")

histogram(vcat(read_dist(10,6)...), label = "", dpi = 250, normalize = true)
title!(L"$n_H = 6$")
xlabel!(L"\lambda")
xlims!(0, 2.5)
ylabel!("p")
savefig("figs/more_samples/6.png")

dists = read_dist(12, 12)

using RandomQuantum
using LinearAlgebra
using StatsBase
using Random

dists_GUE = []
r = RandomQuantum.GUE(2^10)
for i in ProgressBar(1:1024)
	O = sqrt(2)*rand(r)/2^5
	ptranspose!(O, 2^5)
	D = svdvals(O)
	push!(dists_GUE, D)
end

level_spacings(dist) = (d=sort(dist);[(d[i]-d[i-1])/(d[i-1]-d[i-2]) for i in 3:length(d)])

dists12 = read_dist(12, 12)
spacing_dist_12 = []
for dist12 in dists12
	spacing_dist_12 = vcat(spacing_dist_12, level_spacings(dist12))
end

dists6 = read_dist(6, 1)
spacing_dist_6 = []
for dist6 in dists6
	spacing_dist_6 = vcat(spacing_dist_6, level_spacings(dist6))
end

spacing_dist_GUE = []
for distGUE in dists_GUE
	spacing_dist_GUE = vcat(spacing_dist_GUE, level_spacings(distGUE))
end


spacing_dist_perm = []
for distperm in dists_perm
	spacing_dist_perm = vcat(spacing_dist_perm, level_spacings(distperm))
end

using SpecialFunctions
bins = 

histogram(spacing_dist_6[spacing_dist_6 .<= 10], normalize = true, label = L"N = 6, n_H = 1", linetype = :stephist)
histogram!(spacing_dist_12[spacing_dist_12 .<= 10], normalize = true, label = L"N = 12, n_H = 12", linetype = :stephist)
histogram!(spacing_dist_GUE[spacing_dist_GUE .<= 10], normalize = true, linetype = :stephist, label = "GUE")
histogram!(spacing_dist_perm[spacing_dist_perm .<= 10], normalize = true, linetype = :stephist, label = "permutation")

xax = LinRange(0, 5, 100)
a = 1/3
plot!(xax, x->exp(-a)*(a)^x/gamma(x+1))


xlims!(0,5)
title!("Operator level spacing ratio statistics")
xlabel!(L"r")
ylabel!(L"p(r)")

savefig("figs/level_spacing_ratios.png")

histogram(vcat(dists6...), normalize = true, label = "N = 6")
histogram!(vcat(dists12...), normalize = true, label = "N = 12")
histogram!(vcat(dists_GUE...), normalize = true, linetype = :stephist, label = "GUE")
xlims!(0,5)
title!("Operator ES histogram")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")

savefig("figs/ES_histogram.png")


N = 8
nH = 8
H = 1/sqrt(2)*[1 1; 1 -1]
dists_Hperm = []
dim = 2^(N÷2)
for i in ProgressBar(1:1024)
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	O = Base.permutecols!!(Base.permutecols!!(Hf, randperm(dim^2))', randperm(dim^2))'
	ptranspose!(O, dim)
	D = svdvals(O)
	push!(dists_Hperm, D)
end

dist_Hperm_2 = vcat(dists_Hperm...)

dists2 = []
for dist in dists_Hperm
	dists2 = vcat(dists2, level_spacings(dist))
end



histogram(vcat(dists_GUE...), normalize = true, label = "GUE (semicircle)", linestyle = :dash, linewidth = 4, color = :grey)
histogram!(vcat(dist_Hperm_0...), normalize = true, label = "N=8, nH = 1", linetype = :stephist, bins = LinRange(0, 3, 27), color = :red, linewidth = 2)
histogram!(vcat(dist_Hperm_1...), normalize = true, label = "N=8, nH = 2", linetype = :stephist, bins = LinRange(0, 3, 27), color = :purple, linewidth = 2)
histogram!(vcat(dist_Hperm_2...), normalize = true, label = "N=8, nH = 3", linetype = :stephist, bins = LinRange(0, 3, 27), color = :blue, linewidth = 2)

title!(L"ES of $\hat O = (P_1HP_2)$")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")

savefig("figs/PHPconvergence_to_semicircle.png")

histogram(spacing_dist_GUE[spacing_dist_GUE .<= 10], normalize = true)
histogram!(dists[dists .<= 10], normalize = true)

histogram(dists2[dists2 .<= 10], normalize = true)


histogram(spacing_dist_GUE[spacing_dist_GUE .<= 25], normalize = true, label = "GUE (WD)", linestyle = :dash, linewidth = 4, color = :grey)
xlims!(0,5)
histogram!(dists0[dists0 .<= 25], normalize = true, label = "N=8, nH = 1", linetype = :stephist, color = :red, linewidth = 2)
histogram!(dists1[dists1 .<= 25], normalize = true, label = "N=8, nH = 2", linetype = :stephist, color = :purple, linewidth = 2)
histogram!(dists2[dists2 .<= 25], normalize = true, label = "N=8, nH = 3", linetype = :stephist, color = :blue, linewidth = 2, alpha = .7)

title!(L"Level spacing ratios of $\hat O = (P_1HP_2)$")
xlabel!(L"r")
ylabel!(L"p(r)")
savefig("figs/PHPconvergence_to_WD.png")
