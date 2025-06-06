"""Spectral properties of locally generated permutations as a function of depth"""

using Pkg; Pkg.activate("ITensorEnv")
using ITensors
using ProgressBars
using Random
using LinearAlgebra

T = Matrix{ComplexF64}

#create random permutation matrix
rand_perm_mat(d) = T(Base.permutecols!!(Matrix(I, d, d), randperm(d)))

#create a brickwork layer of local permutations
function layer(sites)
    Q1 = Vector{ITensor}()
    for i in 1:Int(floor(length(sites)/3))
        U = rand_perm_mat(8)
        push!(Q1, op(U, sites[3*i-2],sites[3*i-1], sites[3*i]))
    end
    Q2 = Vector{ITensor}()
    for i in 1:Int(floor((length(sites)-1)/3))
        U = rand_perm_mat(8)
        push!(Q2, op(U, sites[3*i-1],sites[3*i], sites[3*i+1]))
    end
    Q3 = Vector{ITensor}()
    for i in 1:Int(floor((length(sites)-2)/3))
        U = rand_perm_mat(8)
        push!(Q3, op(U, sites[3*i], sites[3*i+1], sites[3*i+2]))
    end
    return (Q1,Q2,Q3)
end

#conjugate by layer
function apply_layer(O, l)
    O = apply(l[1], O, apply_dag = true)#apply_day conjugates $O$ by l[1]
    O = apply(l[2], O, apply_dag = true)
    O = apply(l[3], O, apply_dag = true)
    O
end

layers = 15 #max layers
sizes = [6,8,10,12] #system sizes
samples = 2^20 #fixed number of singular values at each size

for s in sizes
    for j in ProgressBar(1:samples*2^(-s))
        dists = [[] for l in 1:layers]
        sites = siteinds(2, s)
        O = MPO(
            sites, 
            [
                [T([1 0; 0 1]) for i in 1:(Int(s/2)-1)];
                [T([0 1; 1 0])]; #X on qubit n/2
                [T([1 0; 0 1]) for i in Int(s/2):s]
            ]
        )
        for l in ProgressBar(1:layers)
            O = apply_layer(O, layer(sites))
            orthogonalize!(O, Int(s/2))
            H,Λ,V = svd(O[Int(s/2)], uniqueinds(O[Int(s/2)], O[Int(s/2)+1]), cutoff = 1E-13)
            dists[l] = diag(Λ)
        end
        for (i,d) in enumerate(dists) #write data
            f = open("data/automaton_MPOs/l$(i)s$(s).dat", append = true, create = true)
            println(f, d)
            close(f)
        end
    end
end