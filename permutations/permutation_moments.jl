using StatsBase
using LinearAlgebra
using ProgressBars
using SimplePartitions
using Random

part_idx(part,j) = argmax(in.(j, part))

include("../tools/permutation_tools.jl")
include("../tools/ES_tools.jl")

t =3 #moment

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
        println(part)
    end
end
sum(arr .== 1) #ideal moment

#https://oeis.org/A094149

#comparison with finite-dize automaton
n = 12
d = 2^n
moms = []
X = kron(
    [
        [[1 0; 0 1] for i in 1:Int(n/2-1)];
        [0 1; 1 0];
        [[1 0; 0 1] for i in 1:Int(n/2)]
    ]...
)
samples = 10
t = 2
for i in ProgressBar(1:samples)
    P = rand_perm_mat(d)
    O = P*X*P'
    ptranspose!(O)
    M = O*O'
    push!(moms,1/d*tr(M^t))
end
mean(moms) #finite-size value
std(moms)/sqrt(length(moms))



# Function to check if a partition is non-crossing
function is_non_crossing(partition)
    for (i, A) in enumerate(partition)
        for (j, B) in enumerate(partition)
            if i < j
                for a in A, b in A
                    for c in B, d in B
                        if a < c < b && !(a < d < b)
                            return false
                        end
                    end
                end
            end
        end
    end
    return true
end

function is_even_odd(s, k)
    a = length([p for p in s if p%2 == 0]) == length([p for p in s if p%2 != 0])
    b = all([(p-1)%(k-1)+1 in s for p in s])
    return a && b
end

# Function to generate and print non-crossing partitions
function even_non_crossing_partitions(n)
    elements = collect(1:n)
    i = 0
    for partition in partitions(elements)
        if all(is_even_odd.(partition))
            println(partition)
        end
    end
end

# Example: Print non-crossing partitions of {1, 2, 3, 4}
n = 6
even_non_crossing_partitions(n)
