#Fig. 8: random automaton circuits. 

include("../tools/imports.jl")
#Data produced by MPOs/automaton_MPOs.jl

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

p = plot(yaxis = :log, palette= :seaborn_dark)
for (n,dat,err) in zip(6:2:12, moms, stds)
	good_points = [l for l in 1:15 if dat[l] > err[l]] 
	#throw out datapoints if the resampling error was too large
	plot!(good_points ./ n, dat[good_points], marker = :o, yerr = err[good_points], label = "N=$(n)", linewidth = 2)
end

ylabel!(L"\frac{\langle \lambda^2 \rangle - 3}{3}")
xlabel!("l/N")
savefig("figures/fig8.svg")

f = open("data/automata_vs_bern/automata_spectrum.txt", read = true)
aut = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

f = open("data/automata_vs_bern/bernoulli_spectrum.txt", read = true)
bern = parse.(Float64, split(read(f, String), "\n")[2:end-1])
close(f)

#inset of Fig. 8
bwidth = .03175
blist = 0:bwidth:3
p = plot(size = (500,200), axesfontsize = 12, guidefontsize=12)
histogram!(aut, bins = blist, normalize = true, label = L"N=12, \text{one sample}", linetype = :stephist, palette = seaborn, linewidth = 2)
histogram!(bern, bins = blist, normalize = true, label = "Bernoulli spectrum", linetype = :stephist, palette = seaborn, linewidth = 2)
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")
yticks!(0:2:12)
savefig("figures/fig8_inset.svg")