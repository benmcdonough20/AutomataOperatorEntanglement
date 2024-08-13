using LinearAlgebra
using RandomQuantum
using StridedViews

include("tools.jl")

N = parse(Int, ARGS[1])
nH = parse(Int, ARGS[2])

BLAS.set_num_threads(parse.(Int, ARGS[3]))

M = 2^14
Ns = collect(6:2:10)
d = (length(Ns), sum(Ns .* (Ns .+ 1) .÷ 2), M)

f = open("data/stor_dat.dat", "a")

samples = Int((2.0)^(-N)*M)

N = 12

r = RandomQuantum.GUE(2)
X = kron([rand(r) for i in 1:N]...)
O = rand_PHP(N, nH)
op = O*X*O'
2^6/sqrt(tr(op'op))
tr(op)

ptranspose!(op)
D = svdvals(op)
using Plots
histogram(D * 4.33, normalize = true)
histogram(level_spacings(D))
xlims!(0, 4)

println(f, ",$(N) $(nH)")
for s in 1:samples
	O = rand_PHP(N, nH)
	op = O*X*O'
	ptranspose!(op)
	println.(f, svdvals(op))
end

close(f)
