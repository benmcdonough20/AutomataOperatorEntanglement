using Pkg; Pkg.activate("itensor_env")
using ITensors
using ProgressBars
using Random
using LinearAlgebra

T = Matrix{ComplexF64}

rand_perm_mat(d) = T(Base.permutecols!!(Matrix(I, d, d), randperm(d)))

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

function apply_layer(O, l)
    O = apply(l[1], O, apply_dag = true)
    O = apply(l[2], O, apply_dag = true)
    O = apply(l[3], O, apply_dag = true)
    O
end

layers = 15
sizes = [6]
samples = [2^12]

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
            f = open("/home/ben/Documents/AutonomaProject/MPOs/dat_perm/l$(i)s$(s).dat", append = true, create = true)
            println(f, d)
            close(f)
        end
    end
end