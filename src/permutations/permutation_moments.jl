include("../tools/imports.jl")

part_idx(part,j) = argmax(in.(j, part))

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