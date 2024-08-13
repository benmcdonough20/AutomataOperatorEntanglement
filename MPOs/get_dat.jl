using Pkg; Pkg.activate("itensor_env")
using Plots
using LaTeXStrings
using StatsBase
using LsqFit

include("../tools/ES_tools.jl")
include("../tools/plotting_defaults.jl")

wigner_CDF(x) = 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)
WD_CDF(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2))
B = 0.233378
A = 2*(π-2)/(4-π)
δWD(x) = 1/4*B*(-2*(1+A)+(2+A)/(1+x)+A/(1+x^2)+(2+A)*atan(x))
WD_CDF_adj(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2)) + δWD(x)

function get_semicircle_dkl(dist)
	distmax = 10
	binnum = 5*sum(dist .< distmax)^(1/3) #x2?
	step = distmax / binnum
	P = [wigner_CDF(i*step)-wigner_CDF((i-1)*step) for i in 1:binnum]
	Q = float.([sum((i-1)*step .<= dist .< i*step) for i in 1:binnum])
	P ./= sum(P)
	Q ./= sum(Q)
	s(x,y) = x == 0 ? 0 : x*log(x/y)
	sum(s.(P,Q))
end


function get_WD_dkl(dist)
	distmax = 5
	binnum = 4*sum(dist .< distmax)^(1/3) #x2?
	#binnum = 100
	step = distmax / binnum
	bins = (0:binnum)*step
	P = float.(StatsBase.fit(Histogram, dist, bins).weights) ./ length(dist)
	Q = [WD_CDF_adj(bins[i+1])-WD_CDF_adj(bins[i]) for (i,b) in enumerate(bins[1:end-1])]
	P ./= sum(P)
	Q ./= sum(Q)
	s(x,y) = x == 0 ? 0 : x*log(x/y)
	sum(s.(P,Q))
end

function read_dat(n, l)
	if n in [6,8]
		#f = open("MPOs/dat/l$(l)s$(n).dat", "r")
		f = open("MPOs/dat_perm/l$(l)s$(n).dat", "r")
	else
		f = open("MPOs/dat_pronto_12/l$(l)s$(n).dat", "r")
	end
	r = read(f, String)
	close(f)
	dats = split(r, "\n")[1:end-1]
	arrs = [parse.(Float64, split(d[5:end-1], ",")) for d in dats]
	arrs = [[arr; zeros(2^n-length(arr))] for arr in arrs]
	arrs
end

dist = read_dat(8, 10)


histogram(vcat(dist...), bins = LinRange(0, 4, 100), normalize = true)
read_dat(6, 1)
read_dat(6,10)[5]
length(vcat(read_dat(8, 1)...))
length(vcat(read_dat(10, 10)...))

layers = 15

#the second moment seems like a good way to monitor the convergence!
cols = [:red, :green, :blue, :orange]
p = plot(yaxis = :log)
for (i,n) in enumerate(6:2:12)
	dkls = [get_WD_dkl(vcat(spacings.(read_dat(n, l))...)) for l in 1:layers]
	plot!((1:layers)/n, dkls, label = "N=$(n)", marker = :x, color = cols[i])
end
display(p)
xlabel!(L"l/N")
yticks!(10.0 .^ [0, -1, -2, -3, -4])
ylabel!(L"D_{\text{KL}}")
title!("Level spacing ratio convergence with circuit depth")
savefig("figs/paper/dkl_evo.svg")

dist = vcat(spacings.(read_dat(6, 3))...)
dist = vcat(spacings.(read_dat(6, 4))...)
histogram(dist[dist .< 15])

using RandomQuantum
r = RandomQuantum.GUE(2^8)
dists = []
for i in 1:2^10
	push!(dists, ES2(rand(r)))
end

get_WD_dkl(vcat(spacings.(dists)...))
rand(r)

cols = [:red, :green, :blue, :orange]
p = plot(yaxis = :log)
for (i,n) in enumerate(6:2:12)
	dkls = [get_semicircle_dkl(vcat(read_dat(n, l)...)) for l in 1:layers]
	plot!((1:layers)*n^b, dkls, label = "N=$(n)", marker = :x, color = cols[i])
	#plot!((1:layers), dkls, label = "N=$(n)", marker = :x, color = cols[i])
end
display(p)
xlabel!(L"l/N")
ylabel!(L"D_{\text{KL}}")
title!("Entanglement spectrum convergence with circuit depth")
savefig("figs/paper/dkl_MP_evo.svg")

rngs = [(3,7), (6,9), (5,9), (7,11)]
dkls = [[get_semicircle_dkl(vcat(read_dat(n, l)...)) for l in 1:layers] for n in 6:2:12]
using CurveFit
params = [linear_fit(r[1]:r[2], log.(dkl[r[1]:r[2]])) for (r,dkl) in zip(rngs,dkls)]
xax = LinRange(2, 8, 100)

for p in params
	plot!(xax, x->exp(p[1]+x*p[2]))
end
display(p)

b
a,b = linear_fit((7:11)/12, log.(dkls[7:11]))



scatter(log.(6:2:12), [log(-p[2]) for p in params],  marker = :X, markersize = 10, xtickfontsize = 16, ytickfontsize = 16, label = "", guidefontsize = 14)
xlabel!(L"\log(N)")
ylabel!(L"\log(\gamma_N)")
a,b = linear_fit(log.(6:2:12), [log(-p[2]) for p in params])
plot!(LinRange(log(6), log(12), 100), x->a+b*x, linewidth = 2, linestyle = :dash, color = :black, label = "")
b
savefig("figs/paper/slope_inset.svg")

a


6*(5/12)


dist = vcat(read_dat(6, 3)...)
h = StatsBase.fit(Histogram,dist[dist .<= 5], nbins = 60)
e = h.edges[1].step.hi
h.weights
xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
inds = (xax .> 0) .& (h.weights .> 0)
plot!(xax[inds], h.weights[inds], xaxis = :log, yaxis = :log)

dist = vcat(read_dat(12, 6)...)
h = StatsBase.fit(Histogram,dist[dist .<= 5], nbins = 200)
h = StatsBase.normalize(h; mode = :pdf)
e = h.edges[1].step.hi
h.weights
xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
inds = (xax .> 0) .& (h.weights .> 0) .& (1:length(xax) .< 10)

@.model(x,p) = -x+p[1]
a = LsqFit.curve_fit(model, log10.(xax[inds]), log10.(h.weights[inds]), [0.0]).param[1]

p = plot()
for (i,n) in enumerate(6:2:12)
	dist = vcat(read_dat(n, round(Int, n/2))...)
	h = StatsBase.fit(Histogram,dist[dist .<= 5], nbins = 200)
	h = StatsBase.normalize(h; mode = :pdf)
	e = h.edges[1].step.hi
	h.weights
	xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
	inds = (xax .> 0) .& (h.weights .> 0)
	plot!(log10.(xax[inds]), log10.(h.weights[inds]), color = cols[i], marker = :x, label = "N=$(n)")
end
display(p)
plot!(log10.(LinRange(0.02, 5, 100)), x-> -x+a, linewidth = 2, color = :black, label = "")
title!(L"Level spacing distribution at $l = N/2$")
xlabel!(L"\log(r)")
ylabel!(L"\log(p(r))")
savefig("figs/paper/oneoverflaw.svg")

function replevel(dist)
	#return sum((dist .< .1) .& (dist .!= Inf) .& (dist .!= NaN))
	return sum((dist .>= .1))
end

cols = [:red, :green, :blue, :orange]
p = plot(yaxis = :log)
for (i,n) in enumerate(6:2:10)
	rep = [replevel(vcat(spacings.(read_dat(n, l))...)) for l in 1:layers]
	plot!((1:layers), rep, label = "N=$(n)", marker = :x, color = cols[i])
end
display(p)
xlabel!(L"l")
yticks!(10.0 .^ [0, -1, -2, -3, -4])
ylabel!(L"D_{\text{KL}}")
title!("Level spacing ratio convergence with circuit depth")
savefig("figs/paper/dkl_evo.svg")

dist = vcat(read_dat(10, 12)...)
get_semicircle_dkl(dist)
plot!(LinRange(0, 2, 100), x->1/π * sqrt(4-x^2))
xlims!(0, 2)
histogram(dist[dist .< 4], normalize = true)

dkls = [get_WD_dkl(vcat(spacings.(read_dat(6, l))...)) for l in 1:layers]

dkls = [get_semicircle_dkl(vcat(read_dat(10, l)...)) for l in 1:layers]
plot(dkls)
dkls

dkls = [mean(vcat(spacings.(read_dat(8, l))...)) for l in 1:layers]
dkls
plot(dkls, yaxis = :log)
hline!([7/4])

histogram(vcat(read_dat(8, 3)...), bins = LinRange(0, 3, 119))

#random permutations

f = open("MPOs/X_data_14.dat", "r")

d = read(f, String)
arrs_str = split(d, "\n")
arrs_str = [a[2:end-1] for a in arrs_str]
arrs = split.(arrs_str, ",")
ideal_dist = vcat([parse.(Float64, a) for a in arrs[1:end-1]]...)

2*(2^17)^(1/3)

function get_ideal_dkl(dist)
	binnum = 101
	step = 3 / binnum
	Q = float.([sum((i-1)*step .<= dist .< i*step) for i in 1:binnum])
	P = float.([sum((i-1)*step .<= ideal_dist .< i*step) for i in 1:binnum])
	P ./= sum(P)
	Q ./= sum(Q)
	s(y,x) = (x == 0 || y == 0) ? 0 : x*log(x/y)
	sum(s.(P,Q))
end

function read_dat(n, l)
	f = open("MPOs/dat_perm/l$(l)s$(n).dat", "r")
	r = read(f, String)
	close(f)
	dats = split(r, "\n")[1:end-1]
	if n in [6,8]
		arrs = [parse.(Float64, split(d[5:end-1], ",")) for d in dats]
	else
		arrs = [parse.(Float64, split(d[2:end-1], ",")) for d in dats]
	end
	arrs = [[arr; zeros(2^n-length(arr))] for arr in arrs]
	arrs
end

d = read_dat(6, 10)

get_ideal_dkl(vcat(d...))

histogram(vcat(d...), bins = LinRange(0, 3, 101), normalize = true)
histogram!(ideal_dist, bins = LinRange(0, 3, 48), normalize = true)

function splitarr(arr)
	sl = Int(length(arr)/4)
	[arr[sl*i+1: sl*(i+1)] for i in 0:3]
end

splitarr(d)
getvar(arr) = std([get_ideal_dkl(v) for v in splitarr(arr)])

dkls = [[get_ideal_dkl(vcat(read_dat(n, l)...)) for l in 1:15] for n in 6:2:12]
stds = [[getvar(vcat(read_dat(n, l)...)) for l in 1:15] for n in 6:2:12]

stds

cols = [:blue, :red, :orange, :green]
p = plot(yaxis = :log)
for (n,c,dat, st) in zip(6:2:12, cols, dkls, stds)
	plot!(collect(1:15) ./  n, dat, label = "", color =c)
	#plot!(collect(1:15), dat, color = c, yerr = st, label = "N = $(n)")
end
display(p)
xlims!(0, 1)
xlabel!(L"$N/l$")
savefig("figs/final_presentation/dkl_scaling_aut_inset.svg")

ylabel!(L"D_{\text{KL}}(\text{Aut. ES})")
xlabel!(L"l")
ylims!(3E-4, 5)
title!("Automaton OES behavior with circuit depth")
savefig("figs/final_presentation/dkl_scaling_aut.svg")

function get_second_moment(dist)
	return mean(dist .^4)
end

function get_one_prob(dist)
	d = round.(dist, digits = 5)
	return sum(d .== 1)/length(d)
end

ideal_one_prob = get_one_prob(ideal_dist)

moms = [[abs(get_second_moment(vcat(read_dat(n, l)...))-3) for l in 1:15] for n in 8:2:12]

p = plot(yaxis = :log)
for (n,c,dat) in zip(8:2:12, cols, moms)
	plot!(collect(1:15) ./ n, dat, color = c)
end
display(p)


oneprobs = [[abs(ideal_one_prob - get_one_prob(vcat(read_dat(n, l)...))) for l in 1:15] for n in 6:2:12]

p = plot(yaxis = :log)
for (n,c,dat) in zip(8:2:12, cols, oneprobs)
	plot!(collect(1:15), dat, color = c, label = "N=$(n)")
end
display(p)
xlabel!(L"l")
ylabel!(L"p(\{1\})")