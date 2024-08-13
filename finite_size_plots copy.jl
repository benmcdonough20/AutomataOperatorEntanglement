using Pkg; Pkg.activate("automata")
using StatsBase
using Plots
using LaTeXStrings
using CurveFit
include("tools/distribution_tools.jl")

RANGE = 7
NUM_BINS = 1000

function read_dist(N, nh)
	if N == 14
		return read_dist_14(nh)
	else
		return read_dist_sub14(N, nh)
	end
end

function read_dist_sub14(N, nh)
	f = open("data/X_N$(N)_nH$(nh).csv", read = true)
	l = read(f, String)
	close(f)
	dat = split(l, ",")
	dists = [parse.(Float64, split(d, "\n")[2:end-1]) for d in dat[2:end-1]]
	overflow = parse.(Float64, split(dat[end], "\n")[2:end-1])
	zeros = parse.(Float64, split(dat[1], "\n")[2:end-1])
	dists, sum(sum(dists)) + length(overflow)
end

function read_dist_14(nh)
	dist1 = read_dat_14(nh, 1)
	dist2 = read_dat_14(nh, 2)
	[dist1[1], dist2[1]], dist1[2] + dist2[2]
end

function read_dat_14(nh, num)
	f = open("data_pronto/data/X_N14_nH$(nh)_$(num).csv", read = true)
	l = read(f, String)
	close(f)
	dat = split(l, ",")
	dist = parse.(Float64, split(dat[2], "\n")[2:end-1])
	overflow = parse.(Float64, split(dat[end], "\n")[2:end-1])
	zeros = parse.(Float64, split(dat[1], "\n")[2:end-1])
	dist, sum(dist) + sum(zeros) + length(overflow)
end

function eCD_dist(dist, binsize, norm = 0)
	if norm == 0
		norm = sum(dist)
	end
	x-> x < binsize ? 0 : sum(dist[1:min(Int(floor(x/binsize)), length(dist))])/norm
end


#note: sample size makes an important difference in this metric, but bin size does not
function get_dist_wass(dist)
	step = RANGE/NUM_BINS*4
	ECDF = eCD_dist(dist, RANGE/NUM_BINS)
	sum([abs(ECDF(i*step) - semicircle_distribution.CDF(i*step)) for i in 1:length(dist)]) * step
end

sc = range(colorant"red", colorant"blue", 5)

p = plot(yaxis = :log, dpi = 250)
for N in 6:2:14
	yerr = 2*[std(get_dist_wass.(read_dist(N, nh)[1])) for nh in 1:N]
	plot!([get_dist_wass(sum(read_dist(N, nh)[1][1:2])) for nh in 1:N], color = sc[(N-4)÷2], yerr = yerr, label = "N=$(N)")
end
title!(L"Wass. dist vs $N$ (X)")
xlabel!(L"n_H")
xticks!(1:14)
ylabel!(L"\Vert P-P_e \Vert_1")
display(p)
savefig("figs/14/wass_dist_X.png")

function rebin(range, num, ecd)
	step = range/num
	return [ecd(step*i)-ecd(step*(i-1)) for i in 1:num]
end

numbins = 1000/10 #results are very different for different number of bins. 100 looks the best
step = RANGE/NUM_BINS

function get_dist_JS(dist)
	numbins = 50
	f = eCD_dist(dist, RANGE/NUM_BINS)
	P = rebin(RANGE, numbins, f)
	Q = rebin(RANGE, numbins, semicircle_distribution.CDF)
	arr = Q .* log.(Q ./ P)
	arr[any.(isnan, arr)] .= 0
	arr[any.(isinf, arr)] .= 0
	return sum(arr)
end

p = plot(yaxis = :log, dpi = 250)

for N in 6:2:14
	yerr = [std(get_dist_JS.(read_dist(N, nh)[1])) for nh in 1:N]
	data = [get_dist_JS(sum(read_dist(N, nh)[1][1:2])) for nh in 1:N]
	yerr = [abs(yerr[i]) > abs(data[i]) ? 0 : yerr[i] for i in 1:length(yerr)]
	plot!([get_dist_JS(sum(read_dist(N, nh)[1][1:2])) for nh in 1:N], color = sc[(N-4)÷2], label = "N=$(N)", yerr = yerr)
end

title!(L"DKL vs $N$ (X)")
xlabel!(L"n_H")
xticks!(1:14)
ylabel!(L"D(P || Q)")
display(p)
savefig("figs/14/dkl_X.png")

function tail_prob(dist)
	binsize = RANGE/NUM_BINS
	num = round(Int,floor(4/binsize))
	return sum(dist[num:end])/sum(dist)
end

p = plot(yaxis = :log, dpi = 250)

for N in 6:2:14
	yerr = [std(tail_prob.(read_dist(N, nh)[1])) for nh in 1:N]
	data = [tail_prob(sum(read_dist(N, nh)[1][1:2])) for nh in 1:N]
	yerr = [abs(yerr[i]) > abs(data[i]) ? 0 : yerr[i] for i in 1:length(yerr)]
	plot!(data, color = sc[(N-4)÷2], label = "N=$(N)", yerr = yerr)
end

title!(L"Tail probability vs $N$ (X)")
xlabel!(L"n_H")
xticks!(1:14)
ylabel!(L"P(\lambda \geq 4)")
display(p)
savefig("figs/14/tail_prob_X.png")