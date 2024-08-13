#Imports
using StatsBase
using LinearAlgebra
using Random
using ProgressBars

#methods
function rand_PHP(N, nH)
    H = exp(1im * π/4 * [0 1; 1 0])
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
    P1 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P2 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P1 * Hf * P2
end

function rand_PHP(N)
    H = 1/sqrt(2) * [1 1; 1 -1]
	#H = exp(1im * π/4 * [0 1; 1 0])
	rand_loc = rand(1:N)
	Hf = insert_op_at_locs(H, N, rand_loc)
    P1 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P2 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P1 * Hf * P2
end

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
	X = insert_op_at_locs(op, N, [loc])
	dim = convert(Int, 2^(N/2))

	#initialize
	for i in ProgressBar(1:samples)
		mats = [rand_PHP(N) for i in 1:nH]
		Mat = product(mats)
		O = Mat*X*Mat'
		ptranspose!(O, dim)
		D = svdvals(O)

		f = open("data/magical_data/$(opname)_N$(N)_nH$(nH)", append = true, create = true)
		println(f, "sample $(i)")
		for d in D
			println(f, d)
		end
		close(f)
	end

end

N = parse(Int,ARGS[1])
nH = parse(Int,ARGS[2])

op = [1 0; 0 -1]
M = 2^14
samples = (2.0)^(-N)*M

run_inst(N, nH, op, round(Int, N/2), samples, "Z")