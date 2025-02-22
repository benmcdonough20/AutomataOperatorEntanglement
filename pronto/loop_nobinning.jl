using LinearAlgebra 
using StatsBase
using Random
using ThreadPinning

BLAS.set_num_threads(parse(Int, ARGS[3]))
println(threadinfo(; blas=true, hints=true))

#methods
function insert_op_at_locs(op, num_qubits, locs)
	return kron([(i in locs) ? op : Matrix(I, 2, 2) for i in 1:num_qubits]...)
end

function ptranspose!(M, hsize)
	conj!(M)
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
N = parse.(Int, ARGS[1])
nH = parse.(Int, ARGS[2])
op = [0.0 1.0; 1.0 0.0]
opname = "X"
H = 1/sqrt(2)*[1.0 1.0; 1.0 -1.0]

#don't modifiy
loc = round(Int, N/2)
X = insert_op_at_locs(op, N, [loc])
dim = convert(Int, 2^(N/2))

M = 2^20
samples = (2.0)^(-N)*M

for i in 1:samples
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	Mat = Base.permutecols!!(Base.permutecols!!(Hf, randperm(dim^2))', randperm(dim^2))'
	O = Mat*X*Mat'
	println("started svd")
	ptranspose!(O, dim)
	D = svdvals(O)
	println("finished, now writing")
	f = open("semicircle_data/$(opname)_N$(N)_nH$(nH).csv",create=true, append=true)	
	println(f, "sample $(i)")
	for d in D
		println(f, d)
	end
	close(f)
end
