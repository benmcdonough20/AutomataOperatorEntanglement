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
using Mmap
using ProgressBars
using StridedViews
using RandomQuantum
using LaTeXStrings
using Plots
include("tools/ES_tools.jl")
include("tools/QI_tools.jl")
include("tools/pauli_tools.jl")

Ns = collect(6:2:10)
M = 2^10
d = (length(Ns), sum(Ns .* (Ns .+ 1) .÷ 2), M)

sum(Ns .* (Ns .+ 1) .÷ 2)

stor_arr = mmap("tmp/PHP.bin", Vector{Float64}, prod(d))
stor_arr_view = sreshape(StridedView(stor_arr), d)

for (i,N) in ProgressBar(enumerate(Ns))
    for nH in ProgressBar(1:N)
        samples = Int((2.0)^(-N)*M)
        sample = sreshape(stor_arr_view[i,nH,:], (samples,2^N))
        for s in 1:samples
            O = rand_PHP(N, nH)
            sample[s,:] .= ES2(O, 2^(N÷2))
        end
    end
end

f = open("tmp/PHP.bin", "w+")
d = (length(Ns), sum(Ns .* (Ns .+ 1) .÷ 2), M)

stor_arr = mmap("tmp_PHP.bin", Vector{Float64}, prod(d), shared = true)
stor_arr_view = sreshape(StridedView(stor_arr), d)

stor_arr_view[1,1,:] .= randn(2^20)

n = 8
dists1 = []
dists2 = []
dists3 = []
r = RandomQuantum.GUE(2)
X1 = kron([rand(r) for i in 1:8]...)
X1 .*= sqrt(2^n/tr(X1'X1))
#tr(X) #GUE product is approximately TRACELESS
IDmat(N) = Matrix{Float64}(I, N, N)
X2 = kron([IDmat(2^3), rand(r), IDmat(2^4)]...)
X2 -= tr(X2)*IDmat(2^n)/2^n #needs to be traceless to remove large eigenvalues apparently...
X2 .*= 2^4 / sqrt(tr(X2'X2)) #normalization

X3 = kron([IDmat(2^4), [0 1.0; 1.0 0], IDmat(2^3)]...)

rand_perm_mat(n) = Base.permutecols!!(IDmat(2^n), randperm(2^n))

#The skewed level spacing results from the real spectrum! Equivalent to taking the real part of a random unitary
r = RandomQuantum.ClosedHaarEnsemble(2^n)
N = 8
H = 1/sqrt(2)*[1 1; 1 -1]
for i in ProgressBar(1:2048)

    P = rand_perm_mat(n)
    Q1 = rand_perm_mat(n)
    Q2 = rand_perm_mat(n)
    Hf = insert_op_at_locs(H, n, sample(1:N, 5, replace = false))
    O1 = P*Hf*Q1
    O2 = P*Hf*Q2

    M = O1 * X1 * O1'
    push!(dists1, ES2(M))

    M = O1 * X2 * O1'
    push!(dists2, ES2(M))

    M = 1/2*(O1*X3*O1'+1im*O2*X3*O1' -1im*O1*X3*O2'+O2*X3*O2')
    push!(dists3, ES2(M))
end

r = RandomQuantum.ClosedHaarEnsemble(2^n)
dists5 = []
for i in ProgressBar(1:2048)
    U = rand(r)
    ptranspose!(U)
    push!(dists5, svdvals(U))
end

histogram(vcat([level_spacings(d) for d in dists1]...), normalize = true, linetype = :stephist, linewidth = 2, color = :red, label = "random product", dpi = 250)
histogram!(vcat([level_spacings(d) for d in dists2]...), normalize = true, linetype = :stephist, linewidth = 2, color = :blue, label = "random local (traceless, normalized)")
histogram!(vcat([level_spacings(d) for d in dists3]...), normalize = true, linetype = :stephist, linewidth = 2, color = :purple, label = "ficticious ancilla")
histogram!(vcat([level_spacings(d) for d in dists4]...), normalize = true, linetype = :stephist, linewidth = 2, color = :gray, label = "GUE")
histogram!(vcat([level_spacings(d) for d in dists5]...), normalize = true, linetype = :stephist, linewidth = 2, color = :gray, label = "Haar")
xax = LinRange(0,5,1000)
Zb = 8/27
b = 1
plot!(xax, x -> 1 / Zb * (x + x^2)^b / (1 + x + x^2)^(1 + 3 / 2 * b),  linewidth = 2, label = L"$\beta$ = 1", color = :green, linestyle = :dash)
Zb = 4/81 * π/sqrt(3)
b = 2
plot!(xax, x -> 1 / Zb * (x + x^2)^b / (1 + x + x^2)^(1 + 3 / 2 * b),  linewidth = 2, label = L"$\beta$ = 2", color = :orange, linestyle = :dash)

title!("WD spacing ratios")
xlabel!(L"r")
ylabel!(L"p(r)")
xlims!(0,5)
savefig("figs/strategies_for_MP_spacing_ratios.png")

#Like the GOE but has diags with variance sqrt(2) times off-diags
n = 8
dists1 = []
dists2 = []
dists3 = []
dists4 = []
r1 = RandomQuantum.GUE(2^n)
r2 = RandomQuantum.ClosedHaarEnsemble(2^n)
for i in ProgressBar(1:2048)
    O = rand_PHP(n, n÷2)
    push!(dists4, ES2(O))
end

n = 3
M = randn(2^n, 2^n)
for i in 1:2^n
    M[i,i] *= sqrt(2)
end
M


histogram(vcat(dists2...), normalize = true, linetype = :stephist, color = :red, linewidth = 2, dpi = 250, label = L"(PHQ)X(PHQ)^\dagger")
histogram!(vcat(dists4...) .* sqrt(2)/2^4, normalize = true, linetype = :stephist, color = :gray, linewidth = 2, dpi = 250, label = "GUE")

title!(L"ES distribution of $X(\tau)$")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
savefig("figs/Xtau_ES_distribution.png")

histogram(vcat([level_spacings(d) for d in dists1]...), normalize = true, linetype = :stephist, linewidth = 2, color = :red, label = L"real($O$) $\ , \ O \sim$ GUE", dpi = 250)
histogram!(vcat([level_spacings(d) for d in dists2]...), normalize = true, linetype = :stephist, linewidth = 2, color = :blue, label = L"(PHQ)X(PHQ)^\dagger")
histogram!(vcat([level_spacings(d) for d in dists3]...), normalize = true, linetype = :stephist, linewidth = 2, color = :purple, label = L"real($UXU^\dagger$)$ \ , \ U\sim$ Haar")
histogram!(vcat([level_spacings(d) for d in dists4]...), normalize = true, linetype = :stephist, linewidth = 2, color = :gray, label = "GUE")

title!("\"Restricted\" real spacing ratio distribution")
xlabel!(L"r")
ylabel!(L"p(r)")
xlims!(0,5)
savefig("figs/restricted_real_spacing_ratios.png")

plot!(LinRange(0,2, 100), x->1/π * sqrt(4-x^2))
#b = 0.10
#plot!(LinRange(0,10,1000), x->(x+x^2)^b/(1+x+x^2)^(1+3/2*b))
