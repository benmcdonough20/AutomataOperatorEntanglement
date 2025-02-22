#Imports
using StatsBase
using Random
using Kronecker
using MKL

T = Matrix{ComplexF64}

function rand_perm_mat(d::Int)::T
	Base.permutecols!!(Matrix(MKL.I, d, d), randperm(d))
end

function insert_op_at_locs(op::T, num_qubits::Int, locs::Vector{Int})::T
	return kronecker([(i in locs) ? op : Matrix(MKL.I, 2, 2) for i in 1:num_qubits]...)
end

function ptranspose!(M::T, hsize::Int)
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

const N = 14

const opname = "X"

const op::T = [0.0 1.0; 1.0 0.0]
const samples = 256

const X::T = insert_op_at_locs(op, N, [Int(N/2)])
const dim = convert(Int, 2^(N/2))

save_dist = []

#initialize
for i in 1:samples
	Mat = rand_perm_mat(2^N)
	O = Mat*X*Mat'
	ptranspose!(O, dim)
	D = MKL.svdvals(O)
	push!(save_dist, D)
end

f = open("/home/benm2/data/big_perm_data/X_data_14.dat", "w")
for d in save_dist
	println(f, d)
end
close(f)
