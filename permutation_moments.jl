using StatsBase
using LinearAlgebra
using ProgressBars
using SimplePartitions
using Random

part_idx(part,j) = argmax(in.(j, part))

include("tools/permutation_tools.jl")
include("tools/ES_tools.jl")

t = 8
arr = []
for (p1, p2) in ProgressBar(Base.product(all_partitions(t), all_partitions(t)))
    c = [(part_idx(p1, i), part_idx(p2, j)) for (i,j) in mom_chain(t)]
    d = length(p1) + length(p2) - length(Set(c))
    push!(arr, d)
end
sum(arr .== 1)

d = 2^12
moms = []
for i in ProgressBar(1:5)
    #S = rand(d,d) .< 1/d
    P = rand_perm_mat(d)
    ptranspose!(P)
    M = P*P'
    push!(moms,1/d*tr(M^4))
end
mean(moms)