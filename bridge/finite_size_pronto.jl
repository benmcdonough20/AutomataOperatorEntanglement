using LinearAlgebra 
using StatsBase
using Random

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
RANGE = 7.0
NUM_BINS = 1000
N = 10
op = [0.0 1.0; 1.0 0.0]
H = 1/sqrt(2)*[1.0 1.0; 1.0 -1.0]
loc = round(Int, N/2)
X = insert_op_at_locs(op, N, [loc])
dim = convert(Int, 2^(N/2))

#initialize
dist = zeros(Int,NUM_BINS)
num_zeros = 0
overflow = []

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
		push!(overflow, d^2)
	else
		dist[k][Int(floor(d^2*NUM_BINS/RANGE))+1] += 1
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
for o in overflow
	println(f, o)
end
close(f)


for (s,n) in zip(samples,Ns)
    for nH in ProgressBar(1:n)
        run_inst(n, nH, X, round(Int, n/2), s, "X", bootstrap_samples)
    end
end