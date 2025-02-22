using LinearAlgebra
using RandomQuantum

include("tools.jl")

N = parse(Int, ARGS[1])
nH = parse(Int, ARGS[2])

BLAS.set_num_threads(parse.(Int, ARGS[3]))

M = 2^19

f = open("/home/benm2/data/magic_data/M_N$(N)_nH$(nH).dat", "w")

samples = Int((2.0)^(-N)*M)

R = exp(1im * 5π/16 * [0 1; 1 0])
X = R * [1 0 ; 0 -1] * R'
X = kron(X, Matrix{Float64}(I, 2^(N-1), 2^(N-1)))
X .*= sqrt(2^N/tr(X'X))

println(f, "$(N) $(nH)")
for s in 1:samples
	O = rand_PHP(N, nH)
	op = O*X*O'
	ptranspose!(op)
	println.(f, ",")
	println.(f, svdvals(op))
end

close(f)
