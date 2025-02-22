#Imports
using StatsBase
using Random
using ThreadPinning
using Kronecker
using Arrow
using DataFrames
using MKL

#MKL.BLAS.set_num_threads(parse.(Int, ARGS[3]))
MKL.BLAS.set_num_threads(Int(ncputhreads()/Threads.nthreads()))
println(threadinfo(; blas = true, hints = true))

const T = Matrix{ComplexF64}
const Rx = T(exp(1im * π/4 * [0 1; 1 0]))
const H = T(1/sqrt(2) * [1 1; 1 -1])

function rand_perm_mat(d::Int)::T
	Base.permutecols!!(Matrix(MKL.I, d, d), randperm(d))
end

#methods
function rand_PHP(N::Int, nH::Int)::T
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	rand_loc = rand(1:N)
	Rxf = insert_op_at_locs(Rx, N, [rand_loc])
	P1 = rand_perm_mat(2^N)
	P2 = rand_perm_mat(2^N)
    	P1 * Hf * Rxf * P2
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

#Run
#
#tunable parameters
#
#constants
#
const N = parse(Int,ARGS[1])
const nH = parse(Int,ARGS[2])

const opname = "Z"

const op::T = [1.0 0.0; 0.0 -1.0]
const M = 2^20
const samples = Int((2.0)^(-N)*M)

const X::T = insert_op_at_locs(op, N, [Int(N/2)])
const dim = convert(Int, 2^(N/2))

save_dist = DataFrame(zeros(Float64, 2^N, samples), :auto)

#initialize
Threads.@threads for i in 1:samples
	Mat = rand_PHP(N, nH)
	O = Mat*X*Mat'
	ptranspose!(O, dim)
	D = MKL.svdvals(O)
	save_dist[!, i] = D
end

Arrow.write("/home/benm2/data/data_one_Rx/$(opname)_N$(N)_nH$(nH).arrow", save_dist)
