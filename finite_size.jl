#Imports
using Pkg
Pkg.activate("automata/")
using StatsBase
using LinearAlgebra
using Random
using RandomQuantum
using ProgressBars
using Plots
using Colors
using LaTeXStrings
using LsqFit

NUM_MOMENTS = 10

include("tools/distribution_tools.jl")

#methods
function insert_op_at_locs(op, num_qubits, locs)
	return kron([(i in locs) ? op : Matrix(I, 2, 2) for i in 1:num_qubits]...)
end

function ptranspose!(M, hsize)
	conj!(M)
	for i in 0:hsize^2-1
	    for j in 0:hsize^2-1
		in2 = i % hsize
		in1 = i ÷ hsize
		out2 = j % hsize
		out1 = j ÷ hsize
		if in2 < out1
		    tmp = M[i+1,j+1]
		    M[i+1,j+1] = M[in1*hsize+out1+1, in2*hsize + out2+1]
		    M[in1*hsize+out1+1, in2*hsize+out2+1] = tmp
		end
	    end
	end
end

#Run
#
#tunable parameters
function run_inst(N, nH, op, loc, samples, opname, bootstrap_samples)

	#constants
	H = 1/sqrt(2)*[1.0 1.0; 1.0 -1.0]
	X = insert_op_at_locs(op, N, [loc])
	dim = convert(Int, 2^(N/2))

	#initialize
	moments = []

	for k in 1:bootstrap_samples
		moments_ret = zeros(NUM_MOMENTS)
		for i in 1:samples
			rand_locs = sample(1:N, nH, replace = false)
			Hf = insert_op_at_locs(H, N, rand_locs)
			Mat = Base.permutecols!!(Base.permutecols!!(Hf, randperm(dim^2))', randperm(dim^2))'
			O = Mat*X*Mat'
			ptranspose!(O, dim)
			D = svdvals(O)
			for j in 1:NUM_MOMENTS
				moments_ret[j] += mean(D .^ (2j))/samples
			end
		end
		push!(moments, moments_ret)
	end
	moments
end

Ns = [6,8,10]
X = [0.0 1.0 ; 1.0 0.0]
Z = [1.0 0.0; 0.0 -1.0]
bootstrap_samples = 6
samples = [2^6, 2^4, 2^2]

ret_dict = Dict()

for (s,n) in zip(samples,Ns)
	ret_dict[n] = []
    for nH in ProgressBar(1:n)
        push!(ret_dict[n], run_inst(n, nH, Z, round(Int, n/2), s, "X", bootstrap_samples))
    end
end

plot!([ret_dict[10][i][1][3] ./ MP_mom(3) for i in 1:10], yaxis = :log)
plot!([ret_dict[10][i][1][4]./ MP_mom(4) for i in 1:10])
plot!([ret_dict[10][i][1][5]./ MP_mom(5) for i in 1:10])

MP_mom(k) = 1/k*sum([binomial(k, r) * binomial(k, r+1) for r in 0:k-1])