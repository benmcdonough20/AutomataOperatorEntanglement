using LinearAlgebra
using StatsBase
using Random

function ptranspose!(M)
	hsize = Int(sqrt(size(M,1)))
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

function insert_op_at_locs(op, num_qubits, locs)
    return kron([(i in locs) ? op : Matrix(I, 2, 2) for i in 1:num_qubits]...)
end

function rand_PHP(N, nH)
	H = 1/sqrt(2) * [1 1; 1 -1]
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	Base.permutecols!!(Base.permutecols!!(Hf, randperm(2^N))', randperm(2^N))'
end
