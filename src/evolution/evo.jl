include("../tools/imports.jl")

## Make data
d = 2^12
n = Int(log2(d))
samples = 64

rngs = [()->randUnitary(d), ()->rand(GUE(d)), ()->rand(GOE(d))]
norms = [1, sqrt(d), sqrt(d)]

dists_classical = [[] for r in rngs]
for i in ProgressBar(1:samples)
    for (i,(r,nm)) in enumerate(zip(rngs, norms))
        M = convert(Matrix{ComplexF64}, r())
        M .*= 1/nm
        push!(dists_classical[i], ES2(M))
    end
end

f = open("data/classical_vs_evo/unitary.dat", create = true, write = true)
println(f, dists_classical[1])
close(f)

f = open("data/classical_vs_evo/GUE.dat", create = true, write = true)
println(f, dists_classical[2])
close(f)

f = open("data/classical_vs_evo/GOE.dat", create = true, write = true)
println(f, dists_classical[3])
close(f)

X = kron(Matrix(I, Int(2.0^(n/2-1)), Int(2.0^(n/2-1))), [0.0 1.0; 1.0 0.0], Matrix(I, Int(2.0^(n/2)), Int(2.0^(n/2))))

function rand_conj_unitary_herm()
    V = randUnitary(d)
    return V*X*V'
end

function rand_conj_unitary_nonherm()
    V = randUnitary(d)
    return V*U*V'
end

function rand_conj_ortho_real()
    O = randOrthogonal(d)
    return O*X*O'
end

rngs = [rand_conj_unitary_herm, rand_conj_ortho_real]

dists_evo = [[] for r in rngs]
for i in ProgressBar(1:samples)
    for (i,r) in enumerate(rngs)
        M = convert(Matrix{ComplexF64}, r())
        M .*= sqrt(d/tr(M*M'))
        push!(dists_evo[i], ES2(M))
    end
end

f = open("data/classical_vs_evo/hermitian_evo.dat", create = true, write = true)
println(f, dists_evo[1])
close(f)

f = open("data/classical_vs_evo/real_evo.dat", create = true, write = true)
println(f, dists_evo[2])
close(f)

