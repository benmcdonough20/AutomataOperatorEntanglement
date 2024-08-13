using Pkg; Pkg.activate("../automata")
using LinearAlgebra
include("../tools/QI_tools.jl")
include("../tools/ES_tools.jl")
using RandomQuantum
using ProgressBars

T = Matrix{ComplexF64}

function layer(s)
    r = ClosedHaarEnsemble(4)
    Q1 = []
    for i in 1:Int(s/2)
        U = rand(r)
        push!(Q1, U)
    end
    Q2 = []
    for i in 1:Int(s/2)-1
        U = rand(r)
        push!(Q2, U)
    end
    kron(Q1...)*kron([1 0; 0 1],kron(Q2...), [1 0; 0 1])
end

layers = 15
sizes = [10]
samples = [2^8]

for (k,s) in enumerate(sizes)
    for j in ProgressBar(1:samples[k])
        dists = [[] for l in 1:layers]
        O = kron([[[1 0.0im; 0 1] for i in 1:Int(s/2)-1]; [[0 1; 1 0]];[[1 0; 0 1] for i in 1:Int(s/2)]]...)
        for l in ProgressBar(1:layers)
            lay = layer(s)
            R = lay'*O*lay
            Λ = ES2(R)
            dists[l] = Λ
            O .= lay'*O*lay
        end
        for (i,d) in enumerate(dists)
            f = open("/home/ben/Documents/AutonomaProject/MPOs/dat_noItensor/l$(i)s$(s).dat", append = true, create = true)
            println(f, d)
            close(f)
        end
    end
end