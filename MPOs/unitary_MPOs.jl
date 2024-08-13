using Pkg; Pkg.activate("automata")
using LinearAlgebra
include("../tools/QI_tools.jl")
include("../tools/ES_tools.jl")
using ITensors
using RandomQuantum
using Plots
using ProgressBars

T = Matrix{ComplexF64}

function layer(sites)
    r = ClosedHaarEnsemble(4)
    Q1 = Vector{ITensor}()
    for i in 1:Int(length(sites)/2)
        U = rand(r)
        push!(Q1, op(U, sites[2*i-1],sites[2*i]))
    end
    Q2 = Vector{ITensor}()
    for i in 1:Int(length(sites)/2)-1
        U = rand(r)
        push!(Q2, op(U, sites[2*i],sites[2*i+1]))
    end
    return (Q1,Q2)
end

function apply_layer(O, l)
    O = apply(l[1], O, apply_dag = true)
    O = apply(l[2], O, apply_dag = true)
    O
end

layers = 15
sizes = [6, 8, 10, 12]
samples = [2^12, 2^10, 2^8, 2^6]

for (k,s) in enumerate(sizes)
    for j in ProgressBar(1:samples[k])
        dists = [[] for l in 1:layers]
        sites = siteinds(2, s)
        O = MPO(sites, [[T([1 0; 0 1]) for i in 1:(Int(s/2)-1)];[T([0 1; 1 0])];[T([1 0; 0 1]) for i in Int(s/2):s]])
        for l in ProgressBar(1:layers)
            O = apply_layer(O, layer(sites))
            orthogonalize!(O, Int(s/2))
            H,Λ,V = svd(O[Int(s/2)], uniqueinds(O[Int(s/2)], O[Int(s/2)+1]), cutoff = 1E-13)
            dists[l] = diag(Λ)
        end
        for (i,d) in enumerate(dists)
            f = open("/home/ben/Documents/AutonomaProject/MPOs/dat/l$(i)s$(s).dat", append = true, create = true)
            println(f, d)
            close(f)
        end
    end
end