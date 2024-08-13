include("tools/imports.jl")


function pauli_sample(N, ct, method)
    rng = RandomQuantum.GUE(2^N)
    for i in 1:ct
        M = rand(rng)
        for i in 1:4^N
            method(pauli(N, i-1), M, i)
        end
    end
end

#how are the Pauli amplitudes distributed? I suppose they are iid Gaussians
dists = [[],[],[],[]]
pauli_sample(1, 10000, (pauli,M,i)->push!(dists[i], tr(pauli*M)))
print(std.(dists))
N = 4
dists = [[] for i in 1:4^N]
pauli_sample(N, 1000, (pauli,M,i)->push!(dists[i], tr(pauli'*M)))
println(std.(dists))
dist = vcat(dists...)
histogram(real.(dist), normalize = true)
xax = LinRange(-3, 3, 1000)
std(dist)
plot!(xax, x->1/sqrt(2*π*σ^2)*exp(-x^2 / (2*σ^2)), linewidth = 3)

#Calculation check for 2x2 case
d = RandomQuantum.GUE(2)
dists = [[],[],[],[]]
for i in 1:10000
    M = rand(d)
    for (e,d) in zip(reshape(M, 4),dists)
        push!(d, e)
    end
end
std.(real.(dists))