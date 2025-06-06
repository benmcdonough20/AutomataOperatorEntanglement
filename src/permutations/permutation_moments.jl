using StatsBase
using LinearAlgebra
using ProgressBars
using SimplePartitions
using Random

using Pkg;Pkg.activate("automata")

part_idx(part,j) = argmax(in.(j, part))

include("../tools/permutation_tools.jl")
include("../tools/ES_tools.jl")

t = 4 #moment

function tuples_to_part(tup)
    d = Dict()
    for (i,t) in enumerate(tup)
        if !(t in keys(d))
            d[t] = [i]
        else
            d[t] = vcat(d[t], i)
        end
    end
    return [d[k] for k in keys(d)]
end

#partition counting
arr = []
for (p1, p2) in Base.product(all_partitions(t), all_partitions(t))
    c = [(part_idx(p1, i), part_idx(p2, j)) for (i,j) in mom_chain(t)]
    d = length(p1) + length(p2) - length(Set(c))
    push!(arr, d)
    part = tuples_to_part(c)
    if d == 1
        println(p1)
        println(p2)
        println()
    end
end
sum(arr .== 1) #ideal moment

#https://oeis.org/A094149

#comparison with finite-dize automaton
n = 10
d = 2^n
moms = []
X = kron(
    [
        [[1 0; 0 1] for i in 1:Int(n/2-1)];
        [0 1; 1 0];
        [[1 0; 0 1] for i in 1:Int(n/2)]
    ]...
)
samples = 100
t = 2
for i in ProgressBar(1:samples)
    #P = rand_perm_mat(d)
    P = rand(d,d) .< 1/d
    #O = P*X*P'
    ptranspose!(P)
    M = P*P'
    push!(moms,(M^t)[1,1])
end
mean(moms) #finite-size value
std(moms)/sqrt(length(moms))


function is_even_odd(s, k)
    a = length([p for p in s if p%2 == 0]) == length([p for p in s if p%2 != 0])
    return a
end


print([p for p in all_partitions(4)])

# Function to generate all the partitions fulfilling a certain property
function check_property(n)
    i = 0
    for partition in all_partitions(n)
        if all(is_even_odd.(partition, n))
            println(partition)
            i+=1
        end
    end
    println(i)
end

# Example: Print non-crossing partitions of {1, 2, 3, 4}
n = 6
check_property(n)