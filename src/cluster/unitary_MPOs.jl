using Pkg; Pkg.activate("automata")
using ITensors
using ProgressBars
using RandomMatrices
#using ThreadPinning

#BLAS.set_num_threads(Int(ncputhreads()/Threads.nthreads()))

T = Matrix{ComplexF64}

function layer(sites)
    r = Haar(2)
    Q1 = Vector{ITensor}()
    for i in 1:Int(length(sites)/2)
        U = Matrix{ComplexF64}(rand(r, 4))
        push!(Q1, op(U, sites[2*i-1],sites[2*i]))
    end
    Q2 = Vector{ITensor}()
    for i in 1:Int(length(sites)/2)-1
        U = Matrix{ComplexF64}(rand(r, 4))
        push!(Q2, op(U, sites[2*i],sites[2*i+1]))
    end
    return (Q1,Q2)
end

function apply_layer(O, l)
    O = apply(l[1], O, apply_dag = true)
    O = apply(l[2], O, apply_dag = true)
    O
end

const layers = 15
const sizes = [6,8]
const samples = [2^(18-s) for s in sizes]

for (k,s) in enumerate(sizes)
    for j in ProgressBar(1:samples[k])
        dists = []
        sites = siteinds(2, s)
        O = MPO(sites, [[T([1 0; 0 1]) for i in 1:(Int(s/2)-1)];[T([0 1; 1 0])];[T([1 0; 0 1]) for i in Int(s/2):s]])
        for l in 1:layers
            O = apply_layer(O, layer(sites))
            orthogonalize!(O, Int(s/2))
            H,Λ,V = svd(O[Int(s/2)], uniqueinds(O[Int(s/2)], O[Int(s/2)+1]), cutoff = 1E-13)
	    push!(dists, diag(Λ))
        end
        for (i,d) in enumerate(dists)
            f = open("data/unitary_MPOs/l$(i)s$(s).dat", append = true, create = true)
            println(f, d)
            close(f)
        end
    end
end