#Imports
using StatsBase
using Random
using MKL

T = Matrix{ComplexF64}

const N = 14
const dim = convert(Int, 2^(N/2))
samples = 256

save_dist = []

#initialize
for i in 1:samples
	O = rand(2^N, 2^N) .<= 1/2^N
	D = MKL.svdvals(O)
	push!(save_dist, D)
end

f = open("/home/benm2/data/big_perm_data/graph_data_14.dat", "w")
for d in save_dist
	println(f, d)
end
close(f)
