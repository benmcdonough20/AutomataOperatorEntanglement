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

include("tools/ES_tools.jl")

wigner_PDF(x) = 0 <= x <= 2 ? 1/π * sqrt(4-x^2) : 0 #Integration convention: use right sums!
wigner_CDF(x) = 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)

WD_PDF(x) = 27/8 * (x+x^2)/(1+x+x^2)^(5/2)
WD_PDF_adj(x) = WD_PDF(x) + B/(1+x)^2*(1/(x+1/x)-A*1/(x+1/x)^2)

WD_CDF(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2))

B = 0.233378
A = 2*(π-2)/(4-π)

δWD(x) = 1/4*B*(-2*(1+A)+(2+A)/(1+x)+A/(1+x^2)+(2+A)*atan(x))
WD_CDF_adj(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2)) + δWD(x)

function get_data(n, nh)
	Arrow.Table("data_pronto/magical_data/X_N$(n)_nH$(nh).arrow")
end

function get_dist_wass(dist)
	ECDF = ecdf(dist)
	step = 1/(2*length(dist)^(1/3))
	maxstep = ceil(maximum(dist)/step)
	sum([abs(ECDF(i*step) - wigner_CDF(i*step)) for i in 1:maxstep]) * step
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

function get_dist_tail(dist)
	sum(dist .>= 2)/length(dist)
end

function get_dist_crossent(dist)
	abs(1 - 3π^2/16 *mean([wigner_PDF(x) for x in dist]))
end

include("tools/plotting_defaults.jl")

#first plot

sterr(x) = std(x)/sqrt(length(x))
dat = [[get_dist_dkl(vcat(get_data(n, nh)...)) for nh in 1:n] for n in 6:2:12]
stds = [[sterr([get_dist_tail(d) for d in get_data(n, nh)]) for nh in 1:n] for n in 6:2:12]

@. model(x,p) = -x*log(2)+p[1]

c = LsqFit.curve_fit(model, 3:6, log.(dat[end][3:6]), [1.0]).param[1]

a,b = linear_fit(6:2:12, [log(minimum(d)) for d in dat])
#c,d = linear_fit(1:4, log.(dat[end][1:4]))

p1 = plot(yaxis = :log)

#xax = LinRange(1, 10, 100)
#xax = LinRange(1, 10, 100)
xax = LinRange(1, 10, 100)
hline!([exp(a+b * N) for N in 6:2:12], linestyle = :dash, color = seaborn[8], label = L"c_1\exp(-\gamma N)")
plot!(xax, x->exp(c-log(2)*x), linestyle = :dash, color = seaborn[8], label = L"c_2\exp(-n_H\ln(2))")
#plot!(xax, x->exp(c+d*x), linestyle = :dash, color = seaborn[8], label = L"c_2\exp(-\gamma n_H)")

for (i,(d,st)) in enumerate(zip(dat, stds))
	plot!(1:(i+2)*2, d, yaxis = :log, color = seaborn[i], yerr=st, label = "N=$(2*i+4)")
end

title!(L"OES scaling of $Z(\tau)$ to MP in $R_x$-doped automaton circuits")
yticks!(10.0 .^[0,-1,-2,-3, -4, -5])
xticks!(1:12)
xlabel!(L"n_H")
ylabel!(L"D_{\text{KL}}")
display(p1)
savefig("figs/paper/DKLfigs/DKL_Z_Rx.svg")
savefig("/home/ben/Desktop/figs/DKL_X_H_replaceentlabel.svg")

#for the weird one
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

#second plot

dist = []
for i in ProgressBar(1:2^10)
	O = randn(2^10, 2^10)
	vals = svdvals(O)
	histogram(vals)
	spacs = spacings(vals)
	push!(dist, spacs...)
end

function get_dist_dkl(dist)
	distmax = 5
	binnum = 2*sum(dist .< distmax)^(1/3) #x2?
	step = distmax / binnum
	bins = (0:binnum)*step
	P = float.(StatsBase.fit(Histogram, dist, bins).weights) ./ length(dist)
	Q = [WD_CDF_adj(bins[i+1])-WD_CDF_adj(bins[i]) for (i,b) in enumerate(bins[1:end-1])]
	P ./= sum(P)
	Q ./= sum(Q)
	s(x,y) = x == 0 ? 0 : x*log(x/y)
	sum(s.(P,Q))
end

dkl_minerr = get_dist_dkl(dist)

dat = vcat(spacings.(get_data(12, 12))...)
histogram(dat[dat .< 25])
xlims!(0,5)

p = plot()
for n in 6:2:12
	dkls = [get_dist_dkl(vcat(spacings.(get_data(n, nh); prec = 6)...)) for nh in 1:n]
	stderr(x) = std(x)/sqrt(length(x))
	regroup(dist, nb) = [vcat(dist[(i-1)*nb+1:i*nb]...) for i in 1:Int(length(dist)/nb)]
	stds = [stderr(get_dist_dkl.(regroup(vcat(spacings.(get_data(n,nh))), 16))) for nh in 1:n]
	x = 1:n
	#plot!(x .*n^(1/2), dkls, yaxis = :log, label = "", color = seaborn[(n-4)÷2 % length(seaborn)+1], linestyle = :dash)
	#scatter!(x .*n^(1/2), dkls, yaxis = :log, label = "N=$(n)", yerr = stds, color = seaborn[(n-4)÷2 % length(seaborn)+1], marker = :o, markersize = 7)
	plot!(x, dkls, yaxis = :log, label = "", color = seaborn[(n-4)÷2 % length(seaborn)+1], linestyle = :dash)
	scatter!(x, dkls, yaxis = :log, label = "N=$(n)", yerr = stds, color = seaborn[(n-4)÷2 % length(seaborn)+1], marker = :o, markersize = 7)
end
hline!([1E-4], color = :black, linestyle = :dash, label = "GUE (numerics)")
xticks!(1:2:12)
#xlabel!(L"n_H \sqrt{N}")
xlabel!(L"n_H")
ylabel!(L"D_{\text{KL}}(\text{WD})")
#xlims!(2, 20)
title!(L"Level spacing ratio KL divergence of $X(\tau)$")
#savefig("/home/ben/Desktop/figs/spacingratio_KL_divergence.svg")
savefig("/home/ben/Documents/AutonomaProject/figs/paper/level_spacrat_div.svg")
display(p)

dist1 = vcat(spacings.(get_data(8,1); prec = 10)...)
dist2 = vcat(spacings.(get_data(8,2); prec = 10)...)

dist1 = dist1[(dist1 .!= Inf) .& (dist1 .!= NaN)]
dist2 = dist2[(dist2 .!= Inf) .& (dist2 .!= NaN)]

bins = LinRange(0, 2, 40)
w = bins[2]-bins[1]

hist1 = float.([sum((i-1)*w .<= dist1 .< i*w) for i in 1:length(bins)])
hist2 = float.([sum((i-1)*w .<= dist2 .< i*w) for i in 1:length(bins)])

hist1 ./= (w*length(dist1))
hist2 ./= (w*length(dist2))
seaborn
plot(bins .+ w/2, hist1, color = seaborn[1], label = L"N = 8, n_{R_x} = 1", marker = :x, markersize = 7)
plot!(bins .+ w/2, hist2, color = seaborn[4], label = L"N = 8, n_{R_x} = 2", marker = :x, markersize = 7)
plot!(LinRange(0, 2, 100), WD_PDF_adj, label = "WD", color = :black, linestyle = :dash, linewidth = 2)
xlabel!(L"r")
ylabel!(L"p(r)")
savefig("/home/ben/Desktop/figs/WDDKL_inset.svg")

xlims!(0,2)

#oh wow cool, so one H in this circuit is enough for the level spacing statistics
sterr(x) = std(x)/sqrt(length(x))
regroup(dist, nb) = [vcat(dist[(i-1)*nb+1:i*nb]...) for i in 1:Int(length(dist)/nb)]
stds = [sterr(get_dist_dkl.(regroup(vcat(spacings.(get_data(n,1))), 8))) for n in 6:2:12]
dkls = [get_dist_dkl(vcat(spacings.(get_data(n, 1))...)) for n in 6:2:12]

model(x, p) = p[1].*x .+ p[2]

#xdata = log.(collect(6:2:12))
xdata = log.(collect(6:2:12))
ydata = log.(dkls)

fit = LsqFit.curve_fit(model, xdata, ydata, [1.0,1.0])
fit.param

scatter(6:2:12, dkls, yaxis = :log,  label = L"n_H = 1 \ (data)", markersize = 8, xaxis = :log)
plot!(6:2:12, x->exp.(fit.param[1]*log(x)+b), color = :black, label = "fit")

a
b
exp(b)

title!(L"$D_{\text{KL}}$ scaling with system size")
xlabel!(L"N")
ylabel!(L"D_{\text{KL}}")

savefig("figs/paper/DKLatNH1.svg")

xax = exp.(6:2:12)

(a.*xax.+b)

histogram(spacs, normalize = true)
plot!(LinRange(0,5, 100), WD_PDF)


stds = [sterr(mean.(regroup(vcat(spacings.(get_data(n,1))), 32))) for n in 6:2:12]
mean(filter!(x->!isnan(x) & !isinf(x), (vcat(spacings.(get_data(8, 3))...))))

7/4

include("tools/QI_tools.jl")


M = Matrix{Float64}(rand(2^10,2^10).<.5)
M .*= 2^5/sqrt(tr(M'M)-22^2)

vals = svdvals(M)

histogram(vals, normalize = true)
plot!(LinRange(0, 2, 100), x->1/π*sqrt(4-x^2))
2^5