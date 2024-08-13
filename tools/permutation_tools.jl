loc_one(col) = convert(Int,sum([i for i in 1:length(col)] .* col))
rand_perm_mat(N) = Base.permutecols!!(Matrix{Float64}(I, N, N), randperm(N))

function cycle(perm_mat, indices) :: Tuple{Array{Int}, Array{Int}}
    start = first(indices)
    a = start
    n = start
    cycle = []
    while true
        push!(cycle, n)
        n = loc_one(perm_mat[:,a])
        n != start || break
        a = n
    end
    (cycle, indices[(i->!(i in cycle)).(indices)])
end

function cycle_notation(perm_mat)
    cycles = []
    indices = 1:last(size(perm_mat))
    while length(indices) > 0
        new_cycle, new_indices = cycle(perm_mat, indices)
        push!(cycles, new_cycle)
        indices = new_indices
    end
    cycles
end

function cycles_to_perm(cycle_list, N)
    perm = [i for i in 1:N]
    for cyc in cycle_list
        perm = cycle_to_perm(cyc, perm)
    end
    perm
end

function cycle_to_perm(cycle, perm)
    for i in reverse(1:length(cycle)-1)
        perm[cycle[1]],perm[cycle[i+1]] = perm[cycle[i+1]],perm[cycle[1]]
    end
    perm
end

function entanglement_spec(cycle, N)
    subsys_dim = round(Int, sqrt(N))
    U,D,V = tsvd(TensorMap(cycle_to_perm_mat(cycle, N), ℂ^(subsys_dim)*ℂ^(subsys_dim), ℂ^(subsys_dim)*ℂ^(subsys_dim)), (1,3),(2,4))
    diag(D.data).^2
end

cycle_to_perm_mat(cycle, N) = Base.permutecols!!(Matrix{Float64}(I, N, N), cycles_to_perm(cycle, N))