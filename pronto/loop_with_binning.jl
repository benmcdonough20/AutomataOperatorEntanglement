#Imports
using StatsBase
using LinearAlgebra
using Random

#Constants
RANGE = 7
RES = 1000

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
function run_inst(N, nH, op, loc, samples, opname)

	#constants
	H = 1/sqrt(2)*[1.0 1.0; 1.0 -1.0]
	X = insert_op_at_locs(op, N, [loc])
	dim = convert(Int, 2^(N/2))

	#initialize
	dist = zeros(Int,RES)
	num_zeros = 0

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
			else
				dist[round(Int,(d^2/RANGE)*RES)+1] += 1
			end
		end
	end

	f = open("/home/benm2/data/$(opname)_N$(N)_nH$(nH).csv",create=true, write=true)
	println(f, "RANGE ", RANGE, "RES", RES)
	println(f, num_zeros)
	for d in dist
		println(f, d)
	end
	close(f)
end

BLAS.set_num_threads(16)

Ns = [14,16]
X = [0.0 1.0 ; 1.0 0.0]
Z = [1.0 0.0; 0.0 -1.0]
samples = [200, 200, 100, 50, 10, 2]

for (s,n) in zip(samples,Ns)
	for nH in 1:n
		run_inst(n, nH, X, round(Int, n/2), s, "X")
	end
end
