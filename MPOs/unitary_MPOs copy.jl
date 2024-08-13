using Pkg; Pkg.activate("automata")
using LinearAlgebra
include("../tools/QI_tools.jl")
include("../tools/ES_tools.jl")
using ITensors
using RandomQuantum
using Plots
using ProgressBars
using LaTeXStrings

ITensors.disable_warn_order()

function layer(sites)
    r = ClosedHaarEnsemble(4)
    Q1 = Vector{ITensor}()
    for i in 1:Int(length(sites)/2)
        U = rand(r)
        push!(Q1, op(U, s[2*i-1],s[2*i]))
    end
    Q2 = Vector{ITensor}()
    for i in 1:Int(length(sites)/2)-1
        U = rand(r)
        push!(Q2, op(U, s[2*i],s[2*i+1]))
    end
    (Q1,Q2)
end


function apply_layer(O, l)
    O .= product(first(l), O)
    O .= product(last(l), O)
end

layers = 5
samples = 256

sites = 10
s = siteinds(2, sites)
O = op(s,"I")
apply_layer(O, layer(s))
xind = s[1]

X = op(s,"I")
X = product(op([0 1; 1 0], xind), X)

dists = [[] for i in 1:layers]
spacings_dists = [[] for i in 1:layers]
for i in ProgressBar(1:samples)
    O = op(s,"I")
    for l in ProgressBar(1:layers)
        apply_layer(O, layer(s))
        P = product(swapprime(conj(O), 0, 1), product(X, O))
        H,Λ,V = svd(P, s[1:Int(sites/2)]..., (s[1:Int(sites/2)]')...)
        push!(dists[l], diag(Λ)...)
        push!(spacings_dists[l], level_spacings(collect(diag(Λ)))...)
    end
end

#power law test

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
dkls
plot(dkls, yaxis = :log)

Zb = 8/27
b = 1

include("../plotting_defaults.jl")

step = .05
p = plot()
for k in [5:8; 10]
    d = [sum((i-1)*step .<= dists[k] .< i*step) for i in 1:ceil(4/step)]/length(dists[k])
    plot!(d, xaxis = :log, yaxis = :log, label = "depth = $(k)", color = seaborn[k%length(seaborn)+1])
    scatter!(d, xaxis = :log, yaxis = :log, label = "", markersize = 3, marker = :o, color = seaborn[k%length(seaborn)+1], markerstrokewidth = 0)
end
display(p)
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
title!(L"Power-law in DOS around depth $d=7$")
savefig("figs/dos_to_powerlaw.svg")


histogram(1:layers, spacings_dists[1][spacings_dists[1] .< 15], xaxis = :log, yaxis = :log)

dists

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
    plot(p1, p2,p3,p4, layout = (2,2), dpi = 250, plot_title ="ES density and level spacing convergence")
end

gif(anim,  fps = 5)

gif(anim, "figs/7_11_meeting/randunitary_observable_dist_evo.gif", fps = 5)

#threshold tilde r value

#learnability transition??

#bitstring distribution learning over Pauli operators --- learnability of the operator wavefunction?