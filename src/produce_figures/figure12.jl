include("../tools/imports.jl")

N = 2^12
samples = 100
dist = countmap([])

#Let's check that this analysis also works for permutations
X = kron([0 1; 1 0], Matrix(I, Int(N/2), Int(N/2)))
for i in ProgressBar(1:samples)
    M = rand_perm_mat(N)
    #M = M*X*M'
    ptranspose!(M)
    g = SimpleDiGraph(M)
    for v in 1:N
        for cmp in weakly_connected_components(g)
            if v in cmp
                add_countmap!(dist, length(cmp), 1)
            end
        end
    end
end

x = LinRange(1,10, 100)
bar(collect(keys(dist)), collect(values(dist))./sum(values(dist)), label = "emp", dpi = 250, linewidth = 0, color = :blue)
title!("Vertex component sizes in random automata")
psmall = sum([exp(-2k)*(2k)^(k-1)/gamma(k+1) for k in 1:10])
vline!([real(1-psmall)*N], label = L"1-p_{small}")
xlabel!("|C(v)|")
ylabel!("P(|C(v)|)")
savefig("figures/fig12_panel1.svg")

xlims!(1, 10)
ylims!(0,.15)
c = 2
plot!(x, x->exp(-c*x)*(c*x)^(x-1)/gamma(x+1), label = "Poisson branching process", linewidth = 3, )
savefig("figures/fig12_panel2.svg")