include("../tools/imports.jl")

wigner_PDF(x) = 0 <= x <= 2 ? 1/π * sqrt(4-x^2) : 0 #Integration convention: use right sums!
wigner_CDF(x) = 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)

function resample(dist, num_samples)
	sample_size = Int(floor(length(dist)/num_samples))
	return [vcat(dist[(i-1)*sample_size+1:i*sample_size]...) for i in 1:num_samples]
end

function get_data(n, nh)
	#data produced by cluster/sparse_op_csv.jl
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

# Plot Fig. 9
#First panel
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
savefig("figures/fig9_panel1.svg")

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
savefig("figures/fig9_panel1_inset.svg")

#Panel 2: diagonal operator
function get_data(n, nh)
	#data produced by cluser/sparse_op_arrow.jl
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

@.model(x,p) = -x*log(2)+p[1]
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
savefig("figures/fig9_panel2.svg")

#Universal scaling function
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
savefig("figures/fig9_panel2_inset.svg")