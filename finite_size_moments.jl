#Imports
using StatsBase
using LinearAlgebra
using Random
using ProgressBars
using Plots
using Colors
using LaTeXStrings
using LsqFit

RANGE = 7.0
NUM_BINS = 1000

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
	dist = [zeros(Int,NUM_BINS) for i in 1:bootstrap_samples]
	outliers = []
	num_zeros = 0

	for k in 1:bootstrap_samples
		for i in 1:samples
			rand_locs = sample(1:N, nH, replace = false)
			Hf = insert_op_at_locs(H, N, rand_locs)
			Mat = Base.permutecols!!(Base.permutecols!!(Hf, randperm(dim^2))', randperm(dim^2))'
			O = Mat*X*Mat'
			ptranspose!(O, dim)
			D = svdvals(O)
			for d in D
				if d == 0
					num_zeros += 1
				end
				if d^2 >= RANGE
					push!(outliers, d^2)
				else
					dist[k][Int(floor(d^2*NUM_BINS/RANGE))+1] += 1
				end
			end
		end
	end

	f = open("data/$(opname)_N$(N)_nH$(nH).csv",create=true, write=true)
	println(f, "")
	println(f, num_zeros)
	for d in dist
		println(f, ",")
		for v in d
			println(f, v)
		end
	end
		
	println(f, ",")
	for outlier in outliers
		println(f, outlier)
	end
	close(f)
end

Ns = [6,8,10,12]
X = [0.0 1.0 ; 1.0 0.0]
Z = [1.0 0.0; 0.0 -1.0]
bootstrap_samples = 6
samples = [2^8, 2^6, 2^4, 4]

for (s,n) in zip(samples,Ns)
    for nH in ProgressBar(1:n)
        run_inst(n, nH, X, round(Int, n/2), s, "X", bootstrap_samples)
    end
end