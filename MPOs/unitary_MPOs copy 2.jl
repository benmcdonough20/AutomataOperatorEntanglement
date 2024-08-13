using Pkg; Pkg.activate("automata")
using LinearAlgebra
include("../tools/QI_tools.jl")
include("../tools/ES_tools.jl")
include("../tools/permutation_tools.jl")
using ITensors
using RandomQuantum
using Plots
using ProgressBars
using LaTeXStrings

ITensors.disable_warn_order()

X = [0.0 1.0; 1.0 0.0]

TOFF = Matrix{Float64}(I, 8, 8)
TOFF[7:8, 7:8] .= X

CNOT = Matrix{Float64}(I, 4, 4)
CNOT[3:4, 3:4] .= X

gates = [CNOT, TOFF]

function rand_3part(N)
    if N <= 1
        return [0]
    end
    n = rand(2:min(N, 3))
    return vcat(n, rand_3part(N-n))
end

function bilayer(N)
    layer1 = rand_3part(N)
    posts = cumsum(layer1)
    newlst = [0]
    while sum(newlst) < N
        arr = [i for i in 2:min(N-sum(newlst), 3) if !(sum(newlst) + i in posts)]
        if length(arr) == 0
            break
        end
        push!(newlst, rand(arr))
    end
    (layer1[1:end-1], newlst[2:end])
end

function bricklayer(l)
    Q1 = Vector{ITensor}()
    cs = [0;cumsum(l)[1:end-1]]
    for (q,p) in zip(l, cs)
        U = rand_perm_mat(2^q)
        push!(Q1, op(U, [s[p+i+1] for i in 0:q-1]))
    end
    Q1
end

function layer(sites)
    l1, l2 = bilayer(length(sites))
    (bricklayer(l1), bricklayer(l2))
end

perm_ideal_dist = []
for i in ProgressBar(1:10)
    M = rand_perm_mat(2^12)
    push!(perm_ideal_dist, ES2(M)...)
end

function apply_layer(O, l)
    O .= product(first(l), O)
    O .= product(last(l), O)
end

layers = 30
samples = 32

sites = 10
s = siteinds(2, sites)
O = op(s,"I")

dists = [[] for i in 1:layers]
spacings_dists = [[] for i in 1:layers]
for i in ProgressBar(1:samples)
    O = op(s,"I")
    for l in ProgressBar(1:layers)
        apply_layer(O, layer(s))
        H,Λ,V = svd(O, s[1:Int(sites/2)]..., (s[1:Int(sites/2)]')...)
        push!(dists[l], diag(Λ)...)
        push!(spacings_dists[l], level_spacings(collect(diag(Λ)))...)
    end
end

function get_wass_dist(dist)
	binnum = floor(length(dist)^(1/3)*2)*2
	binsize = maximum(dist)/binnum
    C1 = ecdf(real(dist))
    C2 = ecdf(real(perm_ideal_dist))
    sum([abs(C1(x)-C2(x))*binsize for x in binsize*(0:binnum)])
end

wds = [get_wass_dist(d) for d in dists]
wds

histogram(perm_ideal_dist, color = :red, bins = bins, normalize = true, label = "random permutation", yaxis = :log)
ylims!(10^(-2), 15)

bins = LinRange(0,4, 166)
anim = @animate for (l,d,s) in zip(1:layers,dists, spacings_dists)
    p2 = histogram(perm_ideal_dist, color = :red, bins = bins, normalize = true, label = "random permutation", yaxis = :log)
    histogram!(d, color = :blue, normalize = true, label = "brickwork permutations")
    ylims!(10^(-2), 15)
    xlabel!(L"\lambda")
    ylabel!(L"p(\lambda)")
    p1 = plot(wds, yaxis = :log, linewidth = 2, color = :blue, label = "")
    xlabel!("t")
    ylabel!(L"W_1")
    vline!([l], color = :black,label = "")
    plot(p2, p1, layout = (2, 1), plot_title = "Brickwork convergence to random permutation", dpi = 250)
end

gif(anim, fps = 5)
gif(anim, "figs/7_11_meeting/permutation_convergence.gif", fps = 5)
mov(anim, "figs/7_11_meeting/permutation_convergence.mov", fps = 5)

#threshold tilde r value

#learnability transition??

#bitstring distribution learning over Pauli operators --- learnability of the operator wavefunction?