include("../tools/imports.jl")
#data produced by pronto/Haar_MPOs.jl

#Plotting code for Fig. 4

#plotting color palettes
c1 = HSV(seaborn[1])
c2 = HSV(seaborn[4])

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
	plot!((1:layers)/n, dkls, label = "N=$(n)", marker = :o, color = cols[i], markerstrokewidth= 0)
end
title!("Level spacing ratio convergence with circuit depth")
xlabel!(L"l/N")
yticks!(10.0 .^ [0, -1, -2, -3, -4])
ylabel!(L"D_{\text{KL}}")
savefig("figs/paper/dkl_evo.svg")

p = plot(yaxis = :log)
cols = [HSV(c2.h, c2.v, c2.s -.2*i) for i in -2:1]
for (i,n) in enumerate(6:2:12)
	dkls = [get_semicircle_dkl(vcat(read_dat(n, l)...)) for l in 1:layers]
	plot!((1:layers)/n, dkls, label = "N=$(n)", marker = :o, color = cols[i], markerstrokewidth = 0)
end
title!("Entanglement spectrum convergence with circuit depth")
xlabel!(L"l/N")
xlims!(0,2)
ylabel!(L"D_{\text{KL}}")
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

#Unviersal function of l/N? (does not appear in paper figure)
scatter(log.(6:2:12), [log(-p[2]) for p in params],  marker = :X, markersize = 10, xtickfontsize = 16, ytickfontsize = 16, label = "", guidefontsize = 14)
xlabel!(L"\log(N)")
ylabel!(L"\log(\gamma_N)")
a,b = linear_fit(log.(6:2:12), [log(-p[2]) for p in params])
plot!(LinRange(log(6), log(12), 100), x->a+b*x, linewidth = 2, linestyle = :dash, color = :black, label = "")
savefig("figs/paper/slope_inset.svg")

#Scale-invariant power law observed in spacings (Fig. 4 inset)
cols = [HSV(c2.h, c2.v, c2.s -.2*i) for i in -2:1]
p = plot(axesfontsize = 30)
for (i,n) in enumerate(6:2:12)
	dist = vcat(read_dat(n, round(Int, n/2))...)
	h = StatsBase.fit(Histogram,dist, nbins = 100)
	h = StatsBase.normalize(h; mode = :pdf)
	e = h.edges[1].step.hi
	h.weights
	xax = [e*(i+1/2) for i in 0:length(h.edges[1])-2]
	inds = (xax .> 0) .& (h.weights .> 0) .& (1:length(xax) .< 12)
	scatter!(log10.(xax), log10.(h.weights), color = cols[i], label = "", markerstrokewidth = 0)
	@.model(x,p) = -p[2]*x+p[1]
	params = LsqFit.curve_fit(model, log10.(xax[inds]), log10.(h.weights[inds]), [1.0, 0.0]).param
	plot!(log10.(xax), x->-params[2]*x+params[1], color = cols[i], label = "η = $(round(params[2], digits = 2))", linewidth = 2)
end
display(p)
xlabel!(L"$\log(r)$")
ylabel!(L"$p(\log(r))$")
savefig("final_paper_figures/powerlaw_inset.svg")

#Confirm that MP distribution is reached by the final layer
cols = [HSV(c1.h, c1.v, c1.s -.2*i) for i in -2:1]
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