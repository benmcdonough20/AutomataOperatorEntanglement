#Now the thought is the MP statistics diagnoses the convergence to a 4-design,
#or some other way to tell how universal a random gate set is

#What we have found so far is that the ES of a random evolved operator
#converges to MP but only in the case of unitary evolution does the correct level
#spacing appear

#Since PHP is still a universal gate set, maybe this is a better diagnosis.

## First, we want to see / characterize the emergence of MP in P+H and compare to Clifford + T
#in the *evolution* operators

using Pkg
Pkg.activate("automata")
using LinearAlgebra
using ProgressBars
using PyCall
using StridedViews
using Plots
using Random
using StatsBase
using LaTeXStrings
qiskit = pyimport("qiskit.quantum_info")
include("tools/ES_tools.jl")
include("tools/QI_tools.jl")
include("tools/pauli_tools.jl")

Ns = collect(6:2:12)
M = 2^19

sum(Ns .* (Ns .+ 1) .÷ 2)

f = open("data/Cliff.dat", "a")


T = [1 0 ; 0 exp(1im*π/4)]

#Important to note: Evolution operator ES is reached for Clifford gates but only witha lot of T gates

N = 12
rand_state(N) = (ψ = randn(ComplexF64, 2^N); ψ/norm(ψ))

rand_state(2)
using RandomQuantum
r = RandomQuantum.GUE(2)
X = kron([rand(r) for i in 1:N]...)

rand_perm_mat(N) = Base.permutecols!!(Matrix{Float64}(I, N, N), randperm(N))
dists = []
for i in ProgressBar(1:1)
    #cliff1 = qiskit.random_clifford(N).to_operator().data
    #O = cliff1 
    #O = rand_PHP(N, 6)
    O = rand_perm_mat(2^N)
    M = O*X*O'
    ptranspose!(M)
    push!(dists, svdvals(M))
end

#confirmed that the states converge to MP
c = real(sqrt(2^N/tr(X'X)))
histogram((vcat(dists...))*c, normalize = true, linetype = :stephist, dpi = 250, linewidth = 2, color = :blue, label = L"$R(\tau) \ R \sim$ rand. prod.")
xax = LinRange(0,2,1000)
plot!(xax, x-> 1/π * sqrt(4-x^2), color = :red, linewidth = 2, label = "semicircle")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
title!("Entanglement spectrum under Automaton evolution")
savefig("figs/classification/ESunderAut.png")

dist = histogram(vcat([level_spacings(d) for d in dists]...), normalize = true, linetype = :stephist, dpi = 250, linewidth = 2, color = :blue, label = L"$R(\tau) \ R \sim$ rand. prod.")
xlims!(0,5)
Zb = 8/27
b = 1
xax = LinRange(0,5,1000)
plot!(xax, x -> 1 / Zb * (x + x^2)^b / (1 + x + x^2)^(1 + 3 / 2 * b), color = :red, linewidth = 2, label = L"GOE ($\beta$ = 1)")
xlabel!(L"r")
ylabel!(L"p(r)")
title!("Level spacing ratios under Automaton evolution")
savefig("figs/classification/RatiosunderAut.png")

dists = []
n = 12
H = 1/sqrt(2)*[1 1; 1 -1]
for i in ProgressBar(1:25)
    #O = qiskit.random_clifford(n).to_operator().data
    O = rand_perm_mat(2^n)*insert_op_at_locs(H, n, sample(1:n, 4, replace = false))*rand_perm_mat(2^n)
    ψ = kron([rand_state(1) for i in 1:n]...)
    ψ = O*ψ
    ρ = sreshape(StridedView(ψ*ψ'), (2^(n÷2),2^(n÷2),2^(n÷2),2^(n÷2)))
    M = zeros(ComplexF64, 2^(n÷2), 2^(n÷2))
    M = sum([ρ[:,a,:,a] for a in 1:2^(n÷2)])
    push!(dists, sqrt.(svdvals(M)))
end

histogram(vcat(dists...), dpi = 250, label = "")
xlabel!(L"λ")
ylabel!(L"p(λ)")
title!("State entanglement under automaton +H evo.")
savefig("figs/classification/StateunderautomatonH.png")

histogram(vcat([level_spacings(d) for d in dists]...), normalize = true, label = "state ES ratios")
xax = LinRange(0,5,1000)
#Zb = 4/81 * π/sqrt(3)
Zb = 4/81 * π/sqrt(3)
b = 2
xlims!(0,5)
plot!(xax, x -> 1 / Zb * (x + x^2)^b / (1 + x + x^2)^(1 + 3 / 2 * b), color = :red, linewidth = 2, label = L"GUE ($\beta$ = 2)")
xlabel!(L"r")
ylabel!(L"p(r)")
title!("Level spacing ratios under Automaton +H evolution")
savefig("figs/classification/StateratiosunderautH.png")

