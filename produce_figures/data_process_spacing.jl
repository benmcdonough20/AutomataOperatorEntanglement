using Pkg; Pkg.activate("automata")
using Plots
using StatsBase
using LinearAlgebra
using Arrow
using LaTeXStrings
using ProgressBars
using LsqFit
using CurveFit
using RandomQuantum

include("/home/benm/Documents/repos/AutomataOperatorEntanglement/tools/plotting_defaults.jl")
include("/home/benm/Documents/repos/AutomataOperatorEntanglement/tools/ES_tools.jl")

wigner_PDF(x) = 0 <= x <= 2 ? 1/π * sqrt(4-x^2) : 0 #Integration convention: use right sums!
wigner_CDF(x) = 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)

WD_PDF(x) = 27/8 * (x+x^2)/(1+x+x^2)^(5/2)
WD_PDF_adj(x) = WD_PDF(x) + B/(1+x)^2*(1/(x+1/x)-A*1/(x+1/x)^2)

WD_CDF(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2))

B = 0.233378
A = 2*(π-2)/(4-π)

δWD(x) = 1/4*B*(-2*(1+A)+(2+A)/(1+x)+A/(1+x^2)+(2+A)*atan(x))
WD_CDF_adj(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2)) + δWD(x)

function resample(dist, num_samples)
	sample_size = Int(floor(length(dist)/num_samples))
	return [vcat(dist[(i-1)*sample_size+1:i*sample_size]...) for i in 1:num_samples]
end

function get_data(n, nh)
	f = open("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/sparse_H/X_N$(n)_nH$(nh).csv", read = true)
	d = read(f, String)
	dat = split(d, "sample")[2:end]
	[parse.(Float64, split(d, "\n")[2:end-1]) for d in dat]
end

function get_dist_dkl(dist)
	binnum = 2*length(dist)^(1/3)
	range = max(maximum(dist), 2) + .01
	step = range/binnum
	P = [wigner_CDF(i*step)-wigner_CDF((i-1)*step) for i in 1:binnum]
	Q = [sum((i-1)*step .<= dist .< i*step) for i in 1:binnum]/length(dist)
	s(x,y) = x == 0 ? 0 : x*log(x/y)
	sum(s.(P,Q))
end

#first plot
dat = [
	[
		get_dist_dkl(vcat(get_data(n, nh)...))
	for nh in 1:n]
for n in 6:2:12]

stds = [
	[
	std(
		filter(
			isfinite,
			get_dist_dkl.(
				resample(get_data(n, nh), 16)
				)
			)
		)/sqrt(16)
	for nh in 1:n]
for n in 6:2:12]

@. model(x,p) = -x*log(2)+p[1]
a,b = linear_fit(6:2:12, [log(minimum(d)) for d in dat])
c,d = linear_fit(1:3, log.(dat[end][1:3]))

p1 = plot(yaxis = :log)
xax = LinRange(1, 4, 100)
hline!([exp(a+b * N) for N in 6:2:12], linestyle = :dash, color = seaborn[8], label = L"c_1\exp(-\gamma_1 N)")
plot!(xax, x->exp(c+d*x), linestyle = :dash, color = seaborn[8], label = L"c_2\exp(-\gamma_1 n_H)")

for (i,(d,st)) in enumerate(zip(dat, stds))
	good_dp = d .>= st*1.2
	plot!(
		collect(1:(i+2)*2)[good_dp],
		d[good_dp],
		yaxis = :log,
		color = seaborn[i],
		yerr = st[good_dp],
		label = "N=$(2*i+4)",
		linewidth = 2,
		marker = :o,
	)
end

yticks!(10.0 .^[0,-1,-2,-3, -4, -5])
xticks!(1:12)
xlabel!(L"n_H")
ylabel!(L"D_{\text{KL}}(\text{DOS} || \text{MP})")
display(p1)
savefig("final_paper_figures/doped_automata_X_MP_convergence_2.svg")

#inset universal scaling function
p1 = plot(yaxis = :log)
xax = LinRange(-2.5, 0, 100)
hline!([exp(a)], color = seaborn[8], linewidth = 3, linestyle = :dash)
plot!(xax, x->exp(c+d*x), linestyle = :dash, color = seaborn[8], label = L"c_2\exp(-\gamma_1 n_H)", linewidth = 3)
for (i,(dp,st)) in enumerate(zip(dat, stds))
	N = (i+2)*2
	scatter!(
		collect(1:N).- N*(b/d),
		dp .* exp(-b*N),
		yaxis = :log,
		color = seaborn[i],
		yerr = st,
		label = "N=$(N)",
		marker = :o,
		markersize = 14,
		markerstrokewidth = 0
	)
end
display(p1)
savefig("final_paper_figures/X_MP_DKL_scaling_inset.svg")

#Z operator dependence
function get_data(n, nh)
	[copy(col) for col in Arrow.Table("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/sparse_H/Z_N$(n)_nH$(nh).arrow")]
end

dat = [
	[
		get_dist_dkl(vcat(get_data(n, nh)...))
	for nh in 1:n]
for n in 6:2:12]

stds = [
	[
	std(
		filter(
			isfinite,
			get_dist_dkl.(
				resample(get_data(n, nh), 16)
				)
			)
		)/sqrt(16)
	for nh in 1:n]
for n in 6:2:12]

@. model(x,p) = -x*log(2)+p[1]
c = LsqFit.curve_fit(model, 3:6, log.(dat[end][3:6]), [1.0]).param[1]
a,b = linear_fit(6:2:12, [log(minimum(d)) for d in dat])

p1 = plot(yaxis = :log)
xax = LinRange(1, 10, 100)
hline!([exp(a+b * N) for N in 6:2:12], linestyle = :dash, color = seaborn[8], label = L"c_1\exp(-\gamma N)")
plot!(xax, x->exp(c-log(2)*x), linestyle = :dash, color = seaborn[8], label = L"c_2\exp(-n_H\ln(2))")

for (i,(d,st)) in enumerate(zip(dat, stds))
	good_dp = d .>= st*1.2
	plot!(
		collect(1:(i+2)*2)[good_dp],
		d[good_dp],
		yaxis = :log,
		color = seaborn[i],
		yerr = st[good_dp],
		label = "N=$(2*i+4)",
		linewidth = 2,
		marker = :o,
	)
end

yticks!(10.0 .^[0,-1,-2,-3, -4, -5])
xticks!(1:12)
xlims!(3, 12)
ylims!(4*10^(-4),.1)
xlabel!(L"n_H")
ylabel!(L"D_{\text{KL}}(\text{DOS}||\text{MP})")
display(p1)

savefig("final_paper_figures/doped_automata_Z_MP_convergence.svg")

p1 = plot(yaxis = :log)
xax = LinRange(-6, 1, 100)
hline!([exp(a)], color = seaborn[8], linewidth = 3, linestyle = :dash)
plot!(xax, x->exp(c-log(2)*x), linestyle = :dash, color = seaborn[8], label = L"c_2\exp(-\gamma_1 n_H)", linewidth = 3)
for (i,(dp,st)) in enumerate(zip(dat, stds))
	N = (i+2)*2
	scatter!(
		(collect(1:N).+ N*(b/log(2)))[3:end],
		(dp .* exp(-b*N))[3:end],
		yaxis = :log,
		color = seaborn[i],
		yerr = st,
		label = "N=$(N)",
		marker = :o,
		markersize = 14,
		markerstrokewidth =0
	)
end
display(p1)
savefig("final_paper_figures/Z_MP_DKL_scaling_inset.svg")

#change in scaling when H is replaced with Rx
a,b = linear_fit(6:2:12, [log(d[Int(length(d)/2)]) for d in dat])
c,d = linear_fit(0:4, log.(dat[end][6:10]))
#c,d = linear_fit(3:6, log.(dat[end][3:6]))

p = plot(yaxis = :log)
for (i,(dist,st)) in enumerate(zip(dat, stds))
	#n = (2*i+4)
	#x = (3:n) .- (a+b*n)/d .+ c/d
	n = (2*i+4)
	x = (1:n) .- (a+b*n)/d .+ c/d
	#scatter!(x, dist[3:end]/exp(n*b+a), yerr = st[3:end] ./ exp(n*b+a), yaxis = :log, color = seaborn[i], label = "N=$(n)")
	scatter!(x, dist/exp(n*b+a), yerr = st ./ exp(n*b+a), yaxis = :log, color = seaborn[i], label = "N=$(n)")
end

hline!([1], linestyle = :dash, color = seaborn[8], label = "")
xax = LinRange(-5, 0, 100)
plot!(xax, exp.(c .+d .* (xax .- c/d)), linestyle = :dash, color = seaborn[8], label = "")

title!(L"Universal scaling of KL divergence with $N$ for $O = X$")
#title!(L"Universal scaling of $Z(\tau)$ KL divergence for $n_H \geq N/2$")
xlabel!(L"\widetilde x(N)")
ylabel!(L"\widetilde y(N)")
savefig("figs/paper/DKLfigs/universal_X_H.svg")


pp = plot(p2, p3, layout = (2,1), size = (500, 600))
p = plot(p1, pp, layout = (1,2), size = (1200,500))

savefig("figs/paper/weird_DKL_scaling.svg")

display(p)

a,b = linear_fit(6:2:12, [log(d[Int(length(d)/2)]) for d in dat])
c,d = linear_fit(0:4, log.(dat[end][6:10]))
#c,d = linear_fit(3:6, log.(dat[end][3:6]))

p = plot(yaxis = :log)
for (i,(dist,st)) in enumerate(zip(dat, stds))
	n = (2*i+4)
	#x = (1:Int(n/2)) .- (a+b*n)/d .+ c/d
	x = (0:Int(n/2))
	#x2 = (1:2) .- (a+b*n)/d .+ c/d
	#scatter!(x2, dist[1:2]/exp(n*b+a), yaxis = :log, color = seaborn[i], markersize = 5, alpha = .4, label = "")
#/exp(n*b+a)
	scatter!(x, dist[Int(n/2):end] ./ exp(b*n+a), yaxis = :log, color = seaborn[i], label = "N=$(2*i+4)", markersize = 5)
	#scatter!(x, dist[Int(n/2):end], yaxis = :log, color = seaborn[i], label = "N=$(2*i+4)", markersize = 5)
end

#Level spacing plot

#baseline accuracy since WD is not exact
dist = []
for i in ProgressBar(1:2^10)
	O = randn(2^10, 2^10)
	vals = svdvals(O)
	spacs = spacings(vals)
	push!(dist, spacs...)
end

function get_dist_dkl(dist)
	distmax = 5 #power-law distribution is hard to estimate accurately, need to introduce a cutoff
	binnum = 2*sum(dist .< distmax)^(1/3)
	step = distmax / binnum
	bins = (0:binnum)*step
	P = float.(StatsBase.fit(Histogram, dist, bins).weights) ./ length(dist)
	Q = [WD_CDF_adj(bins[i+1])-WD_CDF_adj(bins[i]) for (i,b) in enumerate(bins[1:end-1])]
	P ./= sum(P)
	Q ./= sum(Q)
	s(x,y) = x == 0 ? 0 : x*log(x/y)
	sum(s.(P,Q))
end

dkl_minerr = get_dist_dkl(dist) #should be aroudn 1E-4

function get_data(n, nh)
	[copy(col) for col in Arrow.Table("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/sparse_Rx/X_N$(n)_nH$(nh).arrow")]
end

dat = vcat(spacings.(get_data(12, 12))...)
histogram(dat[dat .< 25])
xlims!(0,5)

p = plot(palette = :seaborn_deep)
for n in 6:2:12
	dkls = [get_dist_dkl(vcat(spacings.(get_data(n, nh); prec = 6)...)) for nh in 1:n]
	stderr(x) = std(x)/sqrt(length(x))
	regroup(dist, nb) = [vcat(dist[(i-1)*nb+1:i*nb]...) for i in 1:Int(length(dist)/nb)]
	stds = [stderr(get_dist_dkl.(regroup(vcat(spacings.(get_data(n,nh))), 16))) for nh in 1:n]
	x = 1:n
	plot!(x, dkls, yaxis = :log,  marker = :o, linewidth = 2, yerr = stds, label = "N=$(n)")
end
hline!([dkl_minerr], color = :black, linestyle = :dash, label = "GUE (numerics)")
xticks!(1:12)
xlabel!(L"n_H")
ylabel!(L"D_{\text{KL}}(\text{spacings} || \text{WD})")
#savefig("/home/ben/Desktop/figs/spacingratio_KL_divergence.svg")
savefig("final_paper_figures/doped_automata_X_WD_convergence.svg")
display(p)

#spacing comparison inset
function get_data(n, nh)
	f = open("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/sparse_H/X_N$(n)_nH$(nh).csv", read = true)
	d = read(f, String)
	dat = split(d, "sample")[2:end]
	[parse.(Float64, split(d, "\n")[2:end-1]) for d in dat]
end
dat1 = vcat(spacings.(get_data(12,12))...)
dat2 = vcat(spacings.([copy(col) for col in Arrow.Table("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/sparse_Rx/X_N12_nH12.arrow")])...)

h1 = StatsBase.fit(Histogram,dat1[dat1 .<= 25], nbins = 400)
h2 = StatsBase.fit(Histogram,dat2[dat2 .<= 25], nbins = 400)

plot(LinRange(0, 5, 500), WD(1), linewidth = 2, color = c2, label = L"\text{WD},\beta = 1")
plot!(LinRange(0, 5, 500), WD(0), linewidth =2, color=c1, label = L"p_{GOE}")

c1 = seaborn[1]
c2 = seaborn[4]

scatter!(collect(h1.edges[1]) .+ (.5*25/length(h1.weights)), h1.weights ./ length(dat1) * length(h1.weights)/25, color = c1, label = L"sparse $H$", markersize = 8)
scatter!(collect(h2.edges[1]) .+ (.5*25/length(h2.weights)), h2.weights ./ length(dat2) * length(h2.weights)/25, color = c2, label = L"sparse $R_x$", markersize = 8)

xlims!(0, 4)
xlabel!(L"r")
ylabel!(L"p(r)")

savefig("final_paper_figures/X_level_spacing_inset.svg")

#repeat the analysis with Z
function get_data(n, nh)
	[copy(col) for col in Arrow.Table("/home/benm/Documents/repos/AutomataOperatorEntanglement/data/sparse_Rx/Z_N$(n)_nH$(nh).arrow")]
end

dat = vcat(spacings.(get_data(12, 12))...)
histogram(dat[dat .< 25])
xlims!(0,5)

p = plot(palette = :seaborn_deep)
for n in 6:2:12
	dkls = [get_dist_dkl(vcat(spacings.(get_data(n, nh); prec = 6)...)) for nh in 1:n]
	stderr(x) = std(x)/sqrt(length(x))
	regroup(dist, nb) = [vcat(dist[(i-1)*nb+1:i*nb]...) for i in 1:Int(length(dist)/nb)]
	stds = [stderr(get_dist_dkl.(regroup(vcat(spacings.(get_data(n,nh))), 16))) for nh in 1:n]
	x = 1:n
	plot!(x, dkls, yaxis = :log,  marker = :o, linewidth = 2, yerr = stds, label = "N=$(n)")
end
hline!([dkl_minerr], color = :black, linestyle = :dash, label = "GUE (numerics)")
xticks!(1:12)
xlabel!(L"n_H")
ylabel!(L"D_{\text{KL}}(\text{spacings} || \text{WD})")
#savefig("/home/ben/Desktop/figs/spacingratio_KL_divergence.svg")
savefig("final_paper_figures/doped_automata_Z_WD_convergence.svg")
display(p)