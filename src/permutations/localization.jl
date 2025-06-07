#This code is not reflected in the paper; The purpose of this notebook
#is to explore the origin of the singular component of the permutation OES
#in the localization of eigenvectors on subtrees of the graph.
#This is a further interesting connection that could have implications for
#for the study of computational complexity.

include("../tools/imports.jl")

#additional dependencies
using Graphs
using Compose
using GraphPlot
using Colors
using Optim

import Cairo, Fontconfig

gplot_kwargs = (arrowlengthfrac = 0.03, plot_size = (30cm,30cm))

function gen_large_cluster(N)
    M = convert.(Float64, rand(N,N) .< 1/N)
    g = SimpleDiGraph(M)

    cmps = weakly_connected_components(g)
    big_vs = convert.(Int,cmps[length.(cmps) .== maximum(length.(cmps))][1])

    (M, M[big_vs, big_vs])
end

function right_eigenspaces(M)
    vals, vecs = eigen(M'*M)
    vals = round.(vals, digits = 10)
    d = Dict()
    for (e,v) in zip(eachcol(vecs), vals)
        if !haskey(d, v)
            d[v] = e
        else
            d[v] = hcat(d[v], e)
        end
    end
    d
end

#each key in is an eigenspace of M'M. Now the goal is to find a unitary transformation minimizing the IPR
comm(A,B) = A*B-B*A

function AntiSym(offdiags::Vector{Float64}, N)
    M = zeros(Float64, N,N)
    k = 1
    for i in 1:N-1, j in i+1:N
            M[i,j] = offdiags[k]
            M[j,i] = -offdiags[k]
            k = k+1
    end
    M
end

function bmat(i, N)
    ei = zeros(Float64, Int(N*(N-1)/2)); ei[i] = 1
    AntiSym(ei, N)
end

#parameterize SO(n)
function SO(offdiags::Vector{Float64}, N)
    exp(AntiSym(offdiags, N))
end

function separate_subgraphs(M, plotvec)
    g = DiGraph(M')
    rel_nodes = Set((1:length(plotvec))[abs.(plotvec) .> .001*mean(abs.(plotvec)*length(plotvec))])
    new_outnodes = Vector{Int}()
    new_innodes = Vector{Int}()
    for n in collect(rel_nodes)
        for v in outneighbors(g, n)
            push!(new_outnodes, v)
        end
    end
    for n in new_outnodes
        push!(rel_nodes, n) 
    end
    new_rel_nodes =  Set((1:length(plotvec))[abs.(M*plotvec) .> .001*mean(abs.(plotvec)*length(plotvec))])
    for n in new_rel_nodes
        for v in inneighbors(g, n)
            push!(new_innodes, v)
        end
    end
    for n in new_innodes
        push!(rel_nodes, n)
    end

    rel_node_list = collect(rel_nodes)
    gp,_ = induced_subgraph(g, collect(rel_nodes))
    gp,rel_node_list
end

Z,M = gen_large_cluster(100)
M

gplot(DiGraph(Z); gplot_kwargs...)
gplot(DiGraph(M); gplot_kwargs...)

d = right_eigenspaces(M)
d[1]
keys(d)


v = d[1.6601454497]
g, nds = separate_subgraphs(M, v)
gplot(g)
is_cyclic(g)

plot(v)
N = size(v)[2]

v = d[1]

function f(α)
    -sum((v*SO(α, N)).^4)
end

α0 = rand(Int(N*(N-1)/2))
f(α0)
α0 = ret.minimizer
ret = optimize(f, α0, NelderMead(), Optim.Options(iterations = 100000)) #better for large system sizes

resvecs = v*SO(ret.minimizer, N)

heatmap(abs.(resvecs), clims = (0,1))
title!(L"$\lambda=1$ after optimization")
xlabel!(L"\psi")
ylabel!(L"|\langle v | \psi \rangle|")
xticks!(1:12)
savefig("../../figures/localized_eigenvectors.svg")
#Note: kernel does not appear to be localized

IPRs_cyclic = []
IPRs_acyclic = []
for k in ProgressBar(keys(d))
    if k != 0
        v = d[k]
        if size(d[k], 2) > 1
            N = size(v)[2]
            α0 = rand(Int(N*(N-1)/2))
            ret = optimize(f, α0, NelderMead(), Optim.Options(iterations = 500000)) #better for large system sizes
            resvecs = v*SO(ret.minimizer, N)
        else
            resvecs = v
            N = 1
        end
        colsacyc = []
        colscyc = []
        for resv in eachcol(resvecs)
            g, nds = separate_subgraphs(M, resv)
            if !is_cyclic(SimpleGraph(g))
                push!(colsacyc, resv)
            else
                push!(colscyc, resv)
            end
        end
        if length(colsacyc) > 0
            IPR = 0
            for c in colsacyc
                IPR += sum(abs.(c) .>= 1E-3)/length(colsacyc)
            end
            push!(IPRs_acyclic, (k, IPR))
        end 
        if length(colscyc) > 0
            IPR = 0
            for c in colscyc
                IPR += sum(abs.(c) .>= 1E-3)/length(colscyc)
            end
            push!(IPRs_cyclic, (k, IPR))
        end 
    end
end

scatter([p[1] for p in IPRs_acyclic], [p[2] for p in IPRs_acyclic], label = "acyclic support", color = :blue, marker = :X, markersize = 5)
scatter!([p[1] for p in IPRs_cyclic], [p[2] for p in IPRs_cyclic], label = "cyclic support", color = :red, marker = :X, markersize = 5)
title!("Weight of eigenvectors on cyclic and acyclic subgraphs")
xlabel!("λ")
ylabel!(L"weight \ (\# nodes)")

colmap = x->colormap("RdBu", 101)[round(Int,(x+1)/2*100+1)]
gplot(DiGraph(M'), nodefillc = [colmap(v) for v in plotvec]; gplot_kwargs...)