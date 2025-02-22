#data produced by pronto/unitary_MPOs.jl
using Pkg; Pkg.activate("../automata")
using Plots
using LaTeXStrings
using StatsBase
using LsqFit
using CurveFit

include("../tools/plotting_defaults.jl")

#plotting color palettes
c1 = HSV(seaborn[1])
c2 = HSV(seaborn[4])

include("../tools/ES_tools.jl")
include("../tools/plotting_defaults.jl")

wigner_CDF(x) = 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)
WD_CDF(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2))
B = 0.233378
A = 2*(π-2)/(4-π)
δWD(x) = 1/4*B*(-2*(1+A)+(2+A)/(1+x)+A/(1+x^2)+(2+A)*atan(x))
WD_CDF_adj(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2)) + δWD(x)

function get_semicircle_dkl(dist)
	distmax = 5 #truncate here
	binnum = 2*sum(dist .< distmax)^(1/3)
	step = distmax / binnum
	print(maximum(dist))
	P = [wigner_CDF(i*step)-wigner_CDF((i-1)*step) for i in 1:binnum]
	Q = float.([sum((i-1)*step .<= dist .< i*step) for i in 1:binnum])
	P ./= sum(P)
	Q ./= sum(Q)
	s(x,y) = x == 0 ? 0 : x*log(x/y)
	sum(s.(P,Q))
end

function get_WD_dkl(dist)
	distmax = 25
	binnum = 2*sum(dist .< distmax)^(1/3) #x2?
	#binnum = 200
	binnum = 50
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
	f = open("data/unitary_MPOs/l$(l)s$(n).dat", "r")
	r = read(f, String)
	close(f)
	dats = split(r, "\n")[1:end-1]
	arrs = [parse.(Float64, split(d[2:end-1], ",")) for d in dats]
	arrs = [[arr; zeros(2^n-length(arr))] for arr in arrs]
	arrs
end

layers = 15

function spacings(dist)
	divs = [dist[i+1]-dist[i] for i in 1:length(dist)-1]
	[divs[i+1]/divs[i] for i in 1:length(divs)-1]
end

#sample data plot
dist = vcat(spacings.(read_dat(6, 15))...)
sort(dist)
length(dist)
histogram(dist[dist .< 25], normalize = true, linetype = :stephist)
plot!(LinRange(0, 25, 500), WD(1))
xlims!(0, 3)

#for error estimation
function resample(dist, sample_size)
	num_samples = floor(length(dist) ÷ sample_size)
	return [vcat(dist[(i-1)*sample_size+1:i*sample_size]...) for i in 1:num_samples]
end

cols = [HSV(c1.h, c1.v, c1.s -.2*i) for i in -2:1]
p = plot(yaxis = :log)
for (i,n) in enumerate(6:2:12)
	dkls = [get_WD_dkl(vcat(spacings.(read_dat(n, l))...)) for l in 1:layers]
	#stds = [std(get_WD_dkl.(resample(read_dat(n,l), 32)))/sqrt(31) for l in 1:layers]
	#if n == 12
	#	stds = [0 for i in 1:layers]
	#end
	plot!((1:layers)/n, dkls, label = "N=$(n)", marker = :o, color = cols[i], markerstrokewidth= 0)
	#stds are within numerical precision of zero... need more samples
end
display(p)
xlabel!(L"l/N")
yticks!(10.0 .^ [0, -1, -2, -3, -4])
ylabel!(L"D_{\text{KL}}")
#title!("Level spacing ratio convergence with circuit depth")
savefig("figs/paper/dkl_evo.svg")

p = plot(yaxis = :log)
cols = [HSV(c2.h, c2.v, c2.s -.2*i) for i in -2:1]
for (i,n) in enumerate(6:2:12)
	dkls = [get_semicircle_dkl(vcat(read_dat(n, l)...)) for l in 1:layers]
	plot!((1:layers)/n, dkls, label = "N=$(n)", marker = :o, color = cols[i], markerstrokewidth = 0)
end
display(p)
xlabel!(L"l/N")
xlims!(0,2)
ylabel!(L"D_{\text{KL}}")
#title!("Entanglement spectrum convergence with circuit depth")
savefig("final_paper_figures/distribution_evolution_unitary_fixed.svg")

#extract slopes
rngs = [(3,7), (6,9), (5,9), (7,11)]
dkls = [[get_semicircle_dkl(vcat(read_dat(n, l)...)) for l in 1:layers] for n in 6:2:12]
params = [linear_fit(r[1]:r[2], log.(dkl[r[1]:r[2]])) for (r,dkl) in zip(rngs,dkls)]
xax = LinRange(0, 2, 100)
for (p,n) in zip(params, 6:2:12)
	plot!(xax, x->exp(p[1]+x*n*p[2]))
end
display(p)

#extract slope
scatter(log.(6:2:12), [log(-p[2]) for p in params],  marker = :X, markersize = 10, xtickfontsize = 16, ytickfontsize = 16, label = "", guidefontsize = 14)
xlabel!(L"\log(N)")
ylabel!(L"\log(\gamma_N)")
a,b = linear_fit(log.(6:2:12), [log(-p[2]) for p in params])
plot!(LinRange(log(6), log(12), 100), x->a+b*x, linewidth = 2, linestyle = :dash, color = :black, label = "")
savefig("figs/paper/slope_inset.svg")

dist = vcat(read_dat(6, 3)...)
h = StatsBase.fit(Histogram,dist[dist .<= 5], nbins = 60)
e = h.edges[1].step.hi
h.weights
xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
inds = (xax .> 0) .& (h.weights .> 0)
plot!(xax[inds], h.weights[inds], xaxis = :log, yaxis = :log)

dist = vcat(read_dat(12, 6)...)
h = StatsBase.fit(Histogram,dist, nbins = 100)
h = StatsBase.normalize(h; mode = :pdf)
e = h.edges[1].step.hi
h.weights
xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
inds = (xax .> 0) .& (h.weights .> 0) .& (1:length(xax) .< 12)
scatter(log10.(xax), log10.(h.weights))

@.model(x,p) = -p[2]*x+p[1]
params = LsqFit.curve_fit(model, log10.(xax[inds]), log10.(h.weights[inds]), [1.0, 0.0]).param
plot!(log10.(xax), x->-params[2]*x+params[1])

p = plot(axesfontsize = 20)
for (i,n) in enumerate(6:2:12)
	dist = vcat(read_dat(n, 15)...)
	h = StatsBase.fit(Histogram,dist, nbins = 50)
	h = StatsBase.normalize(h; mode = :pdf)
	e = h.edges[1].step.hi
	h.weights
	xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
	inds = (xax .> 0) .& (h.weights .> 0)
	plot!(xax[inds], h.weights[inds], color = cols[i], marker = :o, label = "N=$(n)", markerstrokewidth=0)
end
display(p)
plot!(LinRange(0, 2, 200), x->1/pi*sqrt(4-x^2), linewidth = 4, color = :black, label = "", linestyle = :dash)
ylims!(-2.4, 1.1)
savefig("final_paper_figures/MPLateTime.svg")

cols = [HSV(c2.h, c2.v, c2.s -.2*i) for i in -2:1]
p = plot(axesfontsize = 20)
for n in [5,6,7]
	dist = vcat(read_dat(12, round(Int, n))...)
	h = StatsBase.fit(Histogram,dist, nbins = 100)
	h = StatsBase.normalize(h; mode = :pdf)
	e = h.edges[1].step.hi
	h.weights
	xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
	inds = (xax .> 0) .& (h.weights .> 0)
	plot!(log10.(xax), log10.(h.weights), color = cols[1], marker = :o, label = "N=12", markerstrokewidth=0)
end
plot!(log10.(xax), x->-params[2]*x+params[1], color = :black, label=L"ax^\beta", linewidth = 3, linestyle = :dash)

xlabel!(L"\log(r)")
ylabel!(L"\log(p(r))")
savefig("final_paper_figures/oneoverflaw.svg")

#random permutations
function read_dat(n, l)
	f = open("data/perm_MPOs/l$(l)s$(n).dat", "r")
	r = read(f, String)
	close(f)
	dats = split(r, "\n")[1:end-1]
	if n in [6,8]
		arrs = [parse.(Float64, split(d[5:end-1], ",")) for d in dats] #type not properly cast in these, resulting in additional "Any" string
	else
		arrs = [parse.(Float64, split(d[2:end-1], ",")) for d in dats]
	end
	arrs = [[arr; zeros(2^n-length(arr))] for arr in arrs]
	arrs
end

function get_second_moment(dist)
	return mean(dist .^4)
end

function resample(dist, num_samples)
	sample_size = Int(floor(length(dist)/num_samples))
	return [vcat(dist[(i-1)*sample_size+1:i*sample_size]...) for i in 1:num_samples]
end

moms = [[abs(get_second_moment(vcat(read_dat(n, l)...))-3)/3 for l in 1:15] for n in 6:2:12]
stds = [[std(get_second_moment.(resample(read_dat(n,l), 4)))/(sqrt(4*3)) for l in 1:15] for n in 6:2:12]

for n in 1:4
	for l in 1:15
		if stds[n][l] >= moms[n][l]
			print(n,l)
		end
	end
end

moms[4][9]
stds[4][9]

p = plot(yaxis = :log, palette= :seaborn_dark)
for (n,dat,err) in zip(6:2:12, moms, stds)
	good_points = [l for l in 1:15 if dat[l] > err[l]]
	plot!(good_points ./ n, dat[good_points], marker = :o, yerr = err[good_points], label = "N=$(n)", linewidth = 2)
	#scatter!(good_points ./ n, dat[good_points], marker = nothing, yerr = err)
end
display(p)

ylabel!(L"\frac{\langle \lambda^2 \rangle - 3}{3}")
xlabel!("l/N")
savefig("final_paper_figures/local_permutation_moments.svg")

dist_local = vcat(read_dat(12, 15)...)
f = open("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/N14_sparse_H/X_N14_nH0.csv", "r")
r = read(f, String)
close(f)
dats = split(r, "\n")[1:end-1]
arrs = [parse.(Float64, split(d[2:end-1], ",")) for d in dats]
arrs = [[arr; zeros(2^14-length(arr))] for arr in arrs]

S(arrs[2].^2/sum(arrs[2].^2))

dist_bigperm = vcat(arrs...)

binnum = 66
stop = 4
binsize = stop/binnum
hist_local = [sum([(i-1)*binsize <= d < i*binsize for d in dist_local]) for i in 1:binnum]./length(dist_local)
hist_bigperm = [sum([(i-1)*binsize <= d < i*binsize for d in dist_bigperm]) for i in 1:binnum]./length(dist_bigperm)

hist_local
hist_bigperm

#calculate kernel probability
sum(round.(dist_bigperm, digits = 12) .== 0)/length(dist_bigperm)
#estimated kernel probability:
exp(-exp(-1))+2/exp(1) -1

b_range = range(0, 3, length=87)
histogram(dist_bigperm,  linetype = :stephist, bins = b_range, normalize = true, label = L"random permutation, $N=14$", color = c2, aspect_ratio = .1, linewidth = 2)
hline!([(exp(-exp(-1))+2/exp(1)-1)*87/3])
histogram!(dist_local, linetype = :stephist, bins = b_range, normalize =true, label = L"N=12, l=15", color = c1, linewidth = 2, alpha = .7)
(exp(-exp(-1))+2/exp(1)-1)*87/3
savefig("final_paper_figures/big_perm_local_perm_cmp_inset.svg")