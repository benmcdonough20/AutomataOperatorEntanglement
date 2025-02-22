include("tools/pauli_tools.jl")
using Pkg; Pkg.activate("automata")
using RandomQuantum
using LinearAlgebra
using ProgressBars
using Plots
using StatsBase
using ProgressBars

d = 5
da = 3
db = 2

#X = kron(Matrix{Float64}(I, 2^(d-1), 2^(d-1)), [1 0; 0 0])
X = 1/d*Diagonal(exp.(1im .* randn(2^d)))

X .-= pauli(d,0)*tr(pauli(d,0)*X) #remove identity
X .*= 2^(db)/sqrt(tr(X'*X)) #approx. normalization

#normalize this matrix to have trace norm 1
r = RandomQuantum.ClosedHaarEnsemble(2^d)
ps = [pauli(db, k) for k in 0:4^db-1]

function sample_expec(Pi, Pj)
    U = rand(r)
    XU = conj(U*X*U')
    abs(sum([sum(XU .* kron(Pi,l))*conj(sum(XU .* kron(Pj,l))) for l in ps]))^2
end

M = zeros(4^da, 4^da)
samples = 100
pairs = [(i,j) for i in 1:4^da for j in i+1:4^da]
for (i,j) in ProgressBar(pairs)
    Pi = pauli(da, i-1)
    Pj = pauli(da, j-1)
    M[i, j] = (sqrt(mean([sample_expec(Pi, Pj) for i in 1:samples]))-1/2^db)*sqrt(samples) #var ~1/d_b^2 (4^db in this context)?
end

mean(M[M .!= 0])
heatmap(M)

#Covariance matrix estimation
arr = [(U=rand(r);pauli_basis(U*X*U', d)) for i in ProgressBar(1:100)]
M = [[arr[i][k] for i in 1:length(arr)] for k in 1:4^d]

i = rand(0:4^d-1)
j = rand(0:4^d-1)
cor(M[i], M[j])

d = 6
p = pauli(d, 20)
stds = [(U = rand(r); U[1,:]'*p*U[1,:]) for i in 1:1000]
std(stds)
1/2^d

U = rand(r)
pauli_basis(U*X*U', d)

#When zero is chosen, the variance is spikey, i.e. there are definitive patterns in the disitrbution

1/2^d
1/2^da
1/4^db


4^da*4^da

for i in ProgressBar(0:4^da-1)
    Pi = pauli(da, i)
    for j in 0:4^da-1
        Pj = pauli(da, j)
        ret = 0
        for q in 1:samples
        end
        M[i+1,j+1] = ret / samples
    end
end

for i in 1:4^da
    M[i,i] = 1/sqrt(2^d)
end

heatmap(M)
1/sqrt(2^d)

#da = 3, db = 2, min = 0.35, definitely depends on db (the one summed over)
#da = 2, db = 3, min = 0.175

samples = 100
pairs = rand([(i,j) for i in 0:4^da-1, j in 0:4^da-1], 5)
res = []

for (i,j) in ProgressBar(pairs)
    Pi = pauli(da, i)
    Pj = pauli(da, j)

    r = RandomQuantum.ClosedHaarEnsemble(2^d)
    ret = 0
    for q in 1:samples
        U = rand(r)
        ret += abs(sum([tr(U*X*U'*kron(Pi,pauli(db, k)))*tr(U*X'*U'*kron(Pj,pauli(db, k))) for k in 0:4^db-1]))
    end
    push!(res, (i,j,ret/samples))
end

res

vars = []
for D in 1:100
    samp = []
    push!(vars, std([1/D*randn()^2 for n in 1:1000]))
end

plot(vars, xaxis = :log, yaxis = :log)
plot!(1:100, x->sqrt(2)/x)

sum(abs.(randn(ComplexF64, 10000)/sqrt(10000)).^2)

stds = []
for d in 2 .^ (6:12)
    r = RandomQuantum.ClosedHaarEnsemble(d)
    U = rand(r)
    d = 2^10
    push!(stds, std(real.(U[:, 1])))
end

plot(2 .^ (6:12), stds, xaxis = :log, yaxis = :log)
plot!(2 .^ (6:12), x->1/sqrt(x))

ds = [4,5,6,7,8,9]
std_d = []
for d in ProgressBar(ds)
    #P = pauli(d, 20)*sqrt(2^d)
    P = kron([1 0; 0 -1], Matrix{Float64}(I, 2^(d-1), 2^(d-1)))
    stds = []
    for i in 1:100
        r = RandomQuantum.ClosedHaarEnsemble(2^d)
        U = rand(r)
        lst = [U[:,i]' * P * U[:,i] for i in 1:d]
        push!(stds, std(lst))
    end
    push!(std_d, mean(stds))
end

plot(2 .^ds, std_d, xaxis = :log, yaxis = :log)
plot!(2 .^ds, x->1/sqrt(x))

include("tools/pauli_tools.jl")
n = 4
Z = insert_op_at_locs([1 0; 0 -1], n, [1,2,3,4])

rand_PHP(N) = rand_perm_mat(2^N) * insert_op_at_locs([1 1; 1 -1], N, [1]) * rand_perm_mat(2^N)

M = zeros(2^n, 2^n)
for i in ProgressBar(1:1)
    P = rand_PHP(n)
    for i in 0:2^n-1; for j in 0:2^n-1
        M[i+1,j+1] += real(tr(kron(pauli(Int(n/2),i),pauli(Int(n/2),j)) * P*Z*P'))
        O = kron(pauli(Int(n/2),i),pauli(Int(n/2),j))
    end
end
end

heatmap(M)