using LinearAlgebra 
using StatsBase
using Random
using ThreadPinning

BLAS.set_num_threads(ncputhreads())
println(threadinfo(; blas=true, hints=true))

#methods
function insert_op_at_locs(op, num_qubits, locs)
	return kron([(i in locs) ? op : Matrix(I, 2, 2) for i in 1:num_qubits]...)
end

function ptranspose!(M, hsize)
	for i in 0:hsize^2-1
	    for j in 0:hsize^2-1
			in1, in2 = i ÷ hsize, i % hsize
			out1, out2 = j ÷ hsize, j % hsize
			if in2 < out1
				tmp = M[i+1,j+1]
				M[i+1,j+1] = M[in1*hsize+out1+1, in2*hsize + out2+1]
				M[in1*hsize+out1+1, in2*hsize+out2+1] = tmp
			end
	    end
	end
end

#constants
RANGE = 7.0
NUM_BINS = 1000
N = 14
nH = parse.(Int, ARGS[1])
filenum = parse.(Int, ARGS[2])
op = [0.0 1.0; 1.0 0.0]
opname = "X"
H = 1/sqrt(2)*[1.0 1.0; 1.0 -1.0]

#don't modifiy
loc = round(Int, N/2)
X = insert_op_at_locs(op, N, [loc])
dim = convert(Int, 2^(N/2))

#initialize
dist = zeros(Int,NUM_BINS)
overflow = []

rand_locs = sample(1:N, nH, replace = false)
Hf = insert_op_at_locs(H, N, rand_locs)
Mat = Base.permutecols!!(Base.permutecols!!(Hf, randperm(dim^2))', randperm(dim^2))'
O = Mat*X*Mat'
ptranspose!(O, dim)
D = svdvals(O)

num_zeros = 0
for d in D
	if d == 0
		global num_zeros += 1
	end
	if d^2 >= RANGE
		push!(overflow, d^2)
	else
		dist[Int(floor(d^2*NUM_BINS/RANGE))+1] += 1
	end
end

f = open("data/$(opname)_N$(N)_nH$(nH)_$(filenum).csv",create=true, write=true)
println(f, "")
println(f, num_zeros)
println(f, ",")
for d in dist
	println(f, d)
end
println(f, ",")
for o in overflow
	println(f, o)
end
close(f)
