include("../tools/imports.jl")

#Generate data for Fig. 5
samples = 2^6
ss = 12
Px = [0.0 1.0; 1.0 0.0]
d = 2^ss
X = kron([[[1 0; 0 1] for i in 1:Int(ss/2-1)];[Px];[[1 0; 0 1] for i in 1:Int(ss/2)]]...)

dist_bern = []
for i in ProgressBar(1:samples)
    P = (rand(d,d) .<= 1/d).*exp.(rand([1im,-1im], (d,d)))
    ptranspose!(P)
    push!(dist_bern, svdvals(P)...)
end

dist_aut = []
for i in ProgressBar(1:samples)
    P = rand_perm_mat(d);
    O = P*X*P'
    ptranspose!(O)
    push!(dist_aut, svdvals(O)...)
end

dist_perm = []
for i in ProgressBar(1:samples)
    P = rand_perm_mat(d);
    ptranspose!(P)
    push!(dist_perm, svdvals(P)...)
end

#save data
f = open("data/automata_vs_bern/bernoulli_spectrum.txt", write = true, create = true)
println(f, "samples: 64, system size: 12, matrix type: bernoulli")
for d in dist_bern
println(f, d)
end
close(f)

f = open("data/automata_vs_bern/automata_spectrum.txt", write = true, create = true)
println(f, "samples: 64, system size: 12, matrix type: automata")
for d in dist_aut
println(f, d)
end
close(f)

f = open("data/automata_vs_bern/permutation_spectrum.txt", write = true, create = true)
println(f, "samples: 64, system size: 12, matrix type: permutation")
for d in dist_perm
println(f, d)
end
close(f)

#(Not included in paper) Demonstration that OES is the same for all 
#automaton-evolved off-diagonal operators
D = 12
ops = [kron([0 exp(1im * π/8*i); exp(-1im * π/8*i) 0],Matrix(I, 2^(D-1), 2^(D-1))) for i in 0:4]

dists = [[] for op in ops]
samples = 64
for (j,op) in enumerate(ops)
    for i in ProgressBar(1:samples)
        P = rand_perm_mat(2^D);
        M = P*op*P'
        ptranspose!(M)
        push!(dists[j], svdvals(M)...)
    end
end

f = open("data/off_diagonal/spectrum_vs_angle.dat", write = true)
for d in dists
    println(f, convert.(Float64,d))
end
close(f)

op = kron([1 0; 0 -1],Matrix(I, 2^(D-1), 2^(D-1)))
dist = []
samples = 64
for i in ProgressBar(1:samples)
    P = rand_perm_mat(2^D);
    M = P*op*P'
    ptranspose!(M)
    push!(dist, svdvals(M)...)
end

f = open("data/off_diagonal/diagonal.dat", write = true)
for d in dist
    println(f, d)
end
close(f)

#Compute moments
moms = Dict()
sizes = [6,8,10,12,14]
for s in sizes
    moms[s] = [[],[],[],[]]
end

samples = 2^20 .* 2.0 .^ (-sizes)
for (k,s) in ProgressBar(enumerate(sizes))
    X = kronecker(Matrix(I, 2^(s÷2-1), 2^(s÷2-1)), [0 1; 1 0], Matrix(I, 2^(s÷2), 2^(s÷2)))
    for i in ProgressBar(1:samples[k])
        P = rand_perm_mat(2^s)
        P = P'*X*P
        ptranspose!(P)
        B = P*P'
        M = B
        for j in 1:4
            M = M*B
            push!(moms[s][j],1/2^s*tr(M))
        end
    end
end

f = open("../data/rand_perm_moms.dat", "w")
println(f, moms)
close(f)