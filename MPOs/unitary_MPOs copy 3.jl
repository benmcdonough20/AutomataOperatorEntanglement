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
using StatsBase
using Random


ITensors.disable_warn_order()

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

function apply_layer(O, l)
    O .= product(first(l), O)
    O .= product(last(l), O)
end

layers = 30
samples = 8

sites = 12
s = siteinds(2, sites)
O1 = op(s,"I")
O2 = op(s,"I")
X = op(s, "I")
X = product(op([0 1; 1 0],s[Int(sites/2)]), X)

function rand_hadamard(s, nH)
    s_i = rand(1:length(s))
    h_i = StatsBase.sample(1:length(s), nH, replace = false)
    O = op(s, "I")
    O = product(op([1 0; 0 1im], s[s_i]), O)
    for h in h_i
        O = product(op(1/sqrt(2)*[1 1; 1 -1], s[h]), O)
    end
    O
end

dists = [[] for i in 1:layers]
spacings_dists = [[] for i in 1:layers]
for i in ProgressBar(1:samples)
    O1 = op(s,"I")
    O2 = op(s,"I")
    HS = rand_hadamard(s, 4)
    for l in ProgressBar(1:layers)
        l1 = layer(s)
        l2 = layer(s)
        O1 = product(first(l1), O1)
        O1 = product(last(l1), O1)
        O2 = product(first(l2), O2)
        O2 = product(first(l2), O2)
        P = product(product(O1, HS), O2)
        O = product(swapprime(conj(P), 0, 1), product(X, P))
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

anim = @animate for (l,d,s) in zip(1:layers,dists, spacings_dists)
    try
        histogram(d, color = :blue, normalize = true, label = "brickwork permutations")
    catch
    end
    plot!(LinRange(0,2,100), x->1/π*sqrt(4-x^2), linewidth = 3, color = :red)
    xlabel!(L"\lambda")
    ylabel!(L"p(\lambda)")
    ylims!(0,1)
    xlims!(0,6)
end

function get_dist_DKL(dist, CDF, range)
	binnum = floor(length(dist)^(1/3)*2)
	binsize = maximum([maximum(dist),range])/binnum
	P = [sum(0 .<= dist .- binsize*i .< binsize) for i in 0:binnum-1] ./ length(dist)
	Q = [CDF((i+1)*binsize) - CDF(i*binsize) for i in 0:binnum-1]
	arr = Q .* log.(Q ./ P)
	arr[any.(isnan, arr)] .= 0
	arr[any.(isinf, arr)] .= 0
    if sum(arr) < 0; return NaN; end
	return sum(arr)
end

spacrats = [mean(d) for d in spacings_dists]
dkls = [get_dist_DKL(d, wigner_CDF, 2) for d in dists]

b = 1
Zb = 8/27

anim = @animate for (l,d,s) in zip(1:layers,dists, spacings_dists)
    p1 = plot(label = "data")
    try
        p1 = histogram(d, normalize = true, color = :blue, label = "data")
    catch
    end
    plot!(LinRange(0,2,100), x->1/π * sqrt(4-x^2), label = "Semicircle", color = :red, size = (800,600))
    ylims!(0,1)
    xlims!(0,4)
    xlabel!(L"\sqrt{\lambda}")
    ylabel!(L"p(\sqrt{\lambda})")
    p2 = histogram(s, normalize = true, label = "data", color = :blue)
    plot!(LinRange(0,5,100), r->1/Zb * (r+r^2)^b/(1+r+r^2)^(1+3/2*b), color = :red, label = L"WD, $\beta = 1$")
    ylims!(0,1)
    xlims!(0,5)
    xlabel!(L"r")
    ylabel!(L"p(r)")
    p3 = plot(1:layers, dkls, yaxis = :log, label = L"D_{KL}")
    vline!([l], color = :black, label = "")
    xlabel!("t")
    p4 = plot(1:layers, spacrats, label = L"\langle r \rangle")
    hline!([7/4], color = :red, label = L"\langle r \rangle_{GOE} = 7/4")
    vline!([l], color = :black, label = "")
    xlabel!("t")
    plot(p1, p2,p3,p4, layout = (2,2), dpi = 250, plot_title ="Permutation ES convergence")
end

gif(anim, fps = 5)
#gif(anim, "figs/7_11_meeting/AutomatonH_convergence.gif", fps = 5)
#mov(anim, "figs/7_11_meeting/AutomatonH_convergence.mov", fps = 5)

rs = [[] for i in 6:12]
for i in 6:2:10
    f = open("data/spacrats_$(i).csv", "r")
    l = read(f, String)
    close(f)
    dat = split(l, ",")
    dat[1] = dat[1][2:end]
    dat[end] = dat[end][1:end-2]
    d = parse.(Float64,dat)
    rs[(i-4)÷2] = d
end
using LaTeXStrings
plot(rs[1], linewidth = 2, color = :blue, label = "N = 6", dpi = 250)
title!(L"$\langle r \rangle$ vs depth at different system size")
plot!(rs[2], linewidth = 2, color = :purple, label = "N = 8")
plot!(rs[3], linewidth = 2, color = :red, label = "N = 10")
xlabel!("t")
ylabel!(L"\langle r  \rangle")
savefig("figs/7_11_meeting/r_ratio.png")

rs[1] = spacrats
f = open("data/spacrats_6.csv", "w")
println(f, spacrats)
close(f)
#threshold tilde r value

#learnability transition??

#bitstring distribution learning over Pauli operators --- learnability of the operator wavefunction?
plot(spacrats)
plot!(rs[3])