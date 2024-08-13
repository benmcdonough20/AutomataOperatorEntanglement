using Pkg
Pkg.activate("automata")
using RandomQuantum
using LinearAlgebra
using Plots
using ProgressBars
using LaTeXStrings
using StatsBase

include("tools/distribution_tools.jl")
include("tools/plotting_defaults.jl")
include("tools/ES_tools.jl")

##Level spacing ratios plot
d = 2^10
n = Int(log2(d))
samples = 100
haar = RandomQuantum.ClosedHaarEnsemble(d)
X = kron(Matrix(I, 2^(n÷2-1), 2^(n÷2-1)), [0 1; 1 0], Matrix(I, 2^(n÷2), 2^(n÷2)))

function rand_loc_op()
    U = rand(haar)
    U' * X * U
end

function rand_GOE_op()
    M = randn(d,d)
    M + M'
end

rngs = [()->rand(haar), ()->rand(RandomQuantum.GUE(d)), rand_GOE_op, rand_loc_op]
spacings(dist) = [(dist[i+2]-dist[i+1])/(dist[i+1]-dist[i]) for i in 1:length(dist)-2]

dists = [[] for r in rngs]
for i in ProgressBar(1:samples)
    for (i,r) in enumerate(rngs)
        M = r()
        push!(dists[i], ES2(M))
    end
end

xax = LinRange(0, 2, 100)

Zb = 8/27
b = 1
WD(x) = 1 / Zb * (x + x^2)^b / (1 + x + x^2)^(1 + 3 / 2 * b)

plot(xax, x->27/8*(x+x^2)/(1+x+x^2)^(5/2), color = seaborn[2], label = L"WD, $\beta = 1$")
plot!(xax, x->81/4 * sqrt(3)/π * (x+x^2)^2/(1+x+x^2)^(4), color = seaborn[1], label = L"WD, $\beta = 2$")
plot!(xax, x->3/4*(x+1)/(1+x+x^2)^(3/2), label = L"\text{Surmise for GOE}", color = seaborn[3])

labels = ["Haar", "GUE", "GOE", L"U^\dagger X U"]

for (i,dist) in enumerate(dists)
    x = LinRange(0,5, 60)[1:end-1]
    d = spacings.(dist)
    scatter!(x .+ 2.1/120, bin(vcat(d...), 0, 5, 60), m = :X, color = seaborn[i], label = labels[i], markerstrokewidth =0, markersize = 7)
end

xlabel!(L"r")
ylabel!(L"p(r)")
title!("Entanglement spectrum level spacing ratios")
savefig("figs/paper/fig2.svg")

###Inset
xax = LinRange(0, 2, 100)
plot(xax, x -> 1/π * sqrt(4-x^2), label = "MP", linewidth = 3, linestyle = :dash, color = :black, legendfontsize = 20, xlabelfontsize = 20, ylabelfontsize = 20, color_palette = seaborn)
norms = [1.0, sqrt(d)/sqrt(2), sqrt(d)*sqrt(2), 1.0]
for (i,d) in enumerate(dists)
    dist = vcat(d...) ./ norms[i]
   histogram!(dist, linetype = :stephist, normalize = true, label = "", linewidth = 2) 
end
xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p\left(\sqrt{\lambda}\right)")
savefig("figs/paper/inset2.svg")

Zb = 8/27
b = 1
WD(x) = 1 / Zb * (x + x^2)^b / (1 + x + x^2)^(1 + 3 / 2 * b)
plot!(xax, WD, color = p[1], label = L"WD, $\beta = 1$")

Zb = 4/81 * π / sqrt(3)
b = 2
plot!(xax, WD, color = p[2], label = L"WD, $\beta = 2$")

xlabel!(L"r")
ylabel!(L"p(r)")
title!("Entanglement spectrum level spacing ratios")

histogram(distGUE, normalize = true, label = "emperical PDF", linetype = :stephist, linewidth = 3, dpi = 300, color =:blue)
histogram!(distHaar, normalize = true, label = "emperical PDF", linetype = :stephist, linewidth = 3, dpi = 300, color =:blue)
histogram!(distGOE, normalize = true, label = "emperical PDF", linetype = :stephist, linewidth = 3, dpi = 300, color =:blue)
histogram!(distX, normalize = true, label = "emperical PDF", linetype = :stephist, linewidth = 3, dpi = 300, color =:blue)

savefig("figs/paper/fig2.svg")

xlims!(0,6)
semicircle(x) = 1/(2π*sqrt(x))*sqrt(4-x)
plot!(x, x->semicircle(x), linewidth = 3, label = "Semicircle Distribution", color = :red)

title!("Entanglement Spectrum Distribution of GUE Operator")
xlabel!(L"$\lambda$")
savefig("figs/GUE_entanglement_spectrum.png")

dist_evos = []
for i in ProgressBar(5:7)
    dist_evo = []
    Da = 2^i
    Db = 2^(10-i)

    X = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), [0.0 1.0;1.0 0.0], Matrix(I, Db, Db))
    X .*= Db / sqrt(tr(X'*X)) #normalization

    d = RandomQuantum.ClosedHaarEnsemble(Da*Db)

    for i in ProgressBar(1:128)
        RU = rand(d)
        M = RU*X*RU'
        tens = TensorKit.TensorMap(M, ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
        U, D, V = tsvd(tens, (1,3),(2,4))
        push!(dist_evo, (diag(D.data).^2)...)
    end
    push!(dist_evos, dist_evo)
end

#normalized by operator norm. Conjugation by unitary does not change operator norm,
# and ||X||_op = 1.
p1 = plot()
for i in 0:2
    Da = 2^(7+i)
    Db = 2^(7-i)
    dist = dist_evos[i+1]

    β = Db^2/Da^2
    λp = (1+Db/Da)^2
    λm = (1-Db/Da)^2
    emp_bins = LinRange(λm*.85, λp*1.1, 50) 
    δ = emp_bins[2] - emp_bins[1]
    xax = LinRange(λm, λp, 1000)

    label1 = ""
    label2 = ""
    if i == 0
        label1 = L"ES of $UXU^\dagger \times \sqrt{d_b/d_a}$"
        label2 = "Marcenko-Pastur"
    end
    dist_binned = [sum(emp_bins[i] .<= dist .< emp_bins[i+1]) for i in 1:length(emp_bins)-1] / (length(dist)*δ)
    scatter!(emp_bins[1:end-1] .+ δ/2, dist_binned, m = :X, color = p[i+1], label = label1, markerstrokewidth =0, markersize = 7)
    MP(x) = 1/(2π*x*β) * sqrt((λp-x)*(x-λm))
    plot!(xax, MP, color = p[i+1], linestyle = :dash, linewidth = 2, label = label2)
end
display(p1)
title!(L"Concentration of ES of $UXU^\dagger$ to MP")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")
savefig("figs/paper/fig1.svg")


histogram(dist, normalize = true, linetype = :stephist, label = "GUE", linewidth = 3, color = :blue, dpi = 300, bins =emp_bins)
histogram(dist_evo, normalize = true, linetype = :stephist, label = L"$UXU^\dagger$", linewidth = 3, color = :green)#, bins = emp_bins)
xlims!(0, 4)
x = LinRange(0,4, 100)
plot!(x, semicircle, linewidth = 3, label = "Semicircle distribution", linestyle = :dash, color = :red)
title!("Spectrum of GUE operator and local operator \n evolved under Haar ensemble")
xlabel!(L"$\lambda$")
#savefig("figs/local_op_averaged.png")

dist = []
bins = LinRange(0, 7, 33)

Da = 2^5
Db = 2^5

for i in ProgressBar(1:50)
    RP = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    tens = TensorKit.TensorMap(RP, ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
    U, D, V = tsvd(tens, (1,3), (2,4))
    map(λ -> push!(dist, λ), diag(D.data).^2)
end

histogram(dist, label = "numerics", normalize = true, bins = bins)
title!("Entanglement spectrum of random permutation matrix")
xlabel!(L"$\lambda$")
savefig("figs/spectrum_of_rand_permutation.png")


#What if we try with just a random matrix?
R = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), randn(ComplexF64, 2,2), Matrix(I, Db, Db))
R = (R+R')/2

H = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), 1/sqrt(2)*[1 1;1 -1], Matrix(I, Db, Db))

for i in ProgressBar(1:1)
    RP = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    tens = TensorKit.TensorMap(RP*H*RP', ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
    U, D, V = tsvd(tens, (1,3), (2,4))
    map(λ -> push!(dist, λ), diag(D.data))
end

histogram(dist, label = "numerics", normalize = true)
xlims!(0, 2)
title!("Entanglement spectrum of random local \n operator evolved under automata")
xlabel!(L"$\lambda$")
savefig("figs/spectrum_of_rand_evolved_op.png")
#Check the CNOT thing

CNOT = [1.0 0 0 0; 0 1.0 0 0; 0 0 0 1.0; 0 0 1.0 0]
tens = TensorMap(CNOT, ℂ^2*ℂ^2, ℂ^2*ℂ^2)
U,D,V = tsvd(tens, (1,3), (2,4))
D

#Note something interesting

X = [0 1.0; 1.0 0]
mat = kron(X, X)
tens = TensorMap(mat, ℂ^2*ℂ^2, ℂ^2*ℂ^2)
U,D,V = tsvd(tens, (1,3), (2,4))
D

tens2 = TensorKit.permute(tens, (1,3,2,4))
reshape(tens2.data, (4,4))

P = [1.0 0 0 0; 0 0 0 1.0; 0 1.0 0 0; 0 0 1.0 0]
tens = TensorMap(P, ℂ^2*ℂ^2, ℂ^2*ℂ^2)
U,D,V = tsvd(tens, (1,3), (2,4))
D

P = [0 0 0 1.0; 1.0 0 0 0; 0 1.0 0 0; 0 0 1.0 0]
tens = TensorMap(P, ℂ^2*ℂ^2, ℂ^2*ℂ^2)
U,D,V = tsvd(tens, (1,3), (2,4))
D

#Testing for universality
# I question whether the permutation matrices satisfy a basic property expressed in 
# https://www.math.harvard.edu/media/feier.pdf. This is tested here

dist1 = []
dist2 = []
bins = 0:.21:8

Dbig = 2^6
Dsmall = 2^4
Nsmall = Dbig^4/Dsmall^4

RP = Base.permutecols!!(Matrix{Float64}(I, Dbig^2, Dbig^2), Random.randperm(Dbig^2))
tens = TensorKit.TensorMap(RP, ℂ^Dbig*ℂ^Dbig, ℂ^Dbig*ℂ^Dbig)
U, D, V = tsvd(tens, (1,3), (2,4))
dist1 = diag(D.data).^2

for i in ProgressBar(1:1000)
    RP = Base.permutecols!!(Matrix{Float64}(I, Dsmall^2, Dsmall^2), Random.randperm(Dsmall^2))
    tens = TensorKit.TensorMap(RP, ℂ^Dsmall*ℂ^Dsmall, ℂ^Dsmall*ℂ^Dsmall)
    U, D, V = tsvd(tens, (1,3), (2,4))
    map(λ -> push!(dist2, λ), diag(D.data).^2)
end

histogram(dist1, bins = bins, normalize = true, linetype = :stephist, label = "high-dimension permutation", linewidth  = 3, dpi = 250)
histogram!(dist2, bins = bins, normalize = true, linetype = :stephist, label = "low-dimension permutations", linewidth = 3)
title!("Universality of autonoma ES")
xlabel!(L"$\lambda$")
savefig("figs/permutation_universality.png")

#Entanglement spectrum of shift matrices of different sizes

function shift(dims)
    a,b = dims
    return [ ((x-y-1)%a == 0) ? 1.0 : 0.0 for x in 1:a, y in 1:b]
end

L = 4
perm = [i for i in 1:L^2]
M = Matrix{Float64}(I, L^2, L^2);
M = Base.permutecols!!(M, perm)
tens = TensorMap(M, ℂ^L*ℂ^L,ℂ^L*ℂ^L)
U,D,V = tsvd(tens, (1,3), (2,4))
[d for d in diag(D.data) if d > 0.001]
histogram(diag(D.data), bins = 10)


Dsmall = 2;
perm = [3,2,4,1]
print(perm)
RP = Base.permutecols!!(Matrix{Float64}(I, Dsmall^2, Dsmall^2), perm);
tens = TensorKit.TensorMap(RP, ℂ^Dsmall*ℂ^Dsmall, ℂ^Dsmall*ℂ^Dsmall);
reshape(TensorKit.permute(tens, (1,3,2,4)).data, (Dsmall^2, Dsmall^2))
U,D,V = tsvd(tens, (1,3), (2,4));
diag(D.data)


# Compare Heisenberg picture of evolving local random operator and evolving permutation
distR = []
distX = []

nqubits = 10
Da = 2^5
Db = 2^5
# 
Z = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), [1 0;0 -1], Matrix(I, Db, Db))
X = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), [0 1;1 0], Matrix(I, Db, Db))
Y = Z+0.5im*X
#H = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), 1/sqrt(2)*[1 1;1 -1], Matrix(I, Db, Db))
r = RandomQuantum.GUE(2)
d = RandomQuantum.ClosedHaarEnsemble(2^nqubits)

#G = kron([rand(r) for i in 1:nqubits]...)
M = rand(r) 
M -= Matrix(I, 2, 2) * tr(M)/2
tr(M)
G = kron(M, Matrix{Float64}(I, 2^(nqubits-1), 2^(nqubits-1)))
Y .*= Da/sqrt(tr(Y'Y))

include("tools/ES_tools.jl")
include("tools/QI_tools.jl")
ptranspose!(G)
svdvals(G)

#eigenvalues are big?
#skewed distribution from large eigenvalue!


distR = []
for i in ProgressBar(1:10)
    locs = sample(1:nqubits, 1)
    H = insert_op_at_locs(1/sqrt(2)*[1 1; 1 -1], nqubits, locs)
    RP1 = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    RP2 = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    Mat = (RP2*H*RP1)
    tens = TensorKit.TensorMap(Mat*Y*Mat', ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
    U, D, V = tsvd(tens, (1,3), (2,4))
    push!(distR, diag(D.data)...)
end

xax = LinRange(0, 2, 100)

#so local operator really is different from non-local operator! random unitary is enough...
#local operator is the "hardest" to get

distR

histogram(distR, normalize = true)

plot!(xax, x->1/π * sqrt(4-x^2))
xlims!(0, 5)


include("tools/pauli_tools.jl")
dist = []
p = pauli(6, 0)
for i in 1:10000
    M = randn(2^6, 2^6)
    push!(dist, tr(p*M))
end

std(dist)


for i in ProgressBar(1:1)
    P = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    tens = TensorKit.TensorMap(P*X*P', ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
    U, D, V = tsvd(tens, (1,3), (2,4))
    map(λ -> push!(distR, λ), diag(D.data))
end

distX = []
nqubits = 8
for i in ProgressBar(1:1024)
    locs = sample(1:nqubits, 2)
    H = insert_op_at_locs(1/sqrt(2)*[1 1; 1 -1], nqubits, locs)
    RP1 = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    RP2 = Base.permutecols!!(Matrix{Float64}(I, Da*Db, Da*Db), Random.randperm(Da*Db))
    Mat = (RP2*H*RP1)
    tens = TensorKit.TensorMap(Mat*X*Mat', ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
    U, D, V = tsvd(tens, (1,3), (2,4))
    #map(λ -> push!(distX, λ), diag(D.data))
    push!(distX, diag(D.data))
end

sort(distX) #note that there are some very large singular values for large permutation matrices, but
#they occurr infrequently

#Is RP H RP then essentially a random unitary matrix? Are these two things the same?
x = LinRange(0, 2, 100)
WD = x->1/π*sqrt(4-x^2)

histogram(distR, normalize = true, label = "local Z", linetype = :stephist, linewidth = 3, bins = bins, dpi = 250)
histogram!(distX, normalize = true,label = "local X", linetype = :stephist, linewidth = 3, bins = bins)
plot!(x, WD, label = "Semicircle", linewidth = 3, linestyle = :dash, color = :red)
xlims!(0,3)
title!("Evolution of X and Z")
xlabel!(L"$\sqrt{\lambda}$")
ylabel!("PDF")
savefig("figs/X_Z_operator_evo.png")

#finite-size scaling of the zero mass observed above
ssizes = [2^2, 2^3, 2^4, 2^5, 2^6]
dists = []

bins = LinRange(0, 3, 50)

numH = 2

for s in ssizes
    dist = []

    X = LinearAlgebra.kron(Matrix(I, convert(Integer, s/2), convert(Integer, s/2)), [0 1;1 0], Matrix(I, s, s))
    #H = LinearAlgebra.kron(Matrix(I, convert(Integer, s/2), convert(Integer, s/2)), 1/sqrt(2)*[1 1;1 -1], Matrix(I, s, s))
    nqubits = convert(Int, log2(s)*2)
    locs = sample(1:nqubits, numH, replace = false)
    H = insert_op_at_locs(1/sqrt(2)*[1 1; 1 -1], nqubits, locs)

    println((2^6/s)^3)
    for i in 1:convert(Integer, (2^6/s)^2)
        RP1 = Base.permutecols!!(Matrix{Float64}(I, s^2, s^2), Random.randperm(s^2))
        RP2 = Base.permutecols!!(Matrix{Float64}(I, s^2, s^2), Random.randperm(s^2))
        Mat = (RP2*H*RP1)
        tens = TensorKit.TensorMap(Mat*X*Mat', ℂ^s*ℂ^s, ℂ^s*ℂ^s)
        U, D, V = tsvd(tens, (1,3), (2,4))
        map(λ -> push!(dist, λ), diag(D.data))
    end
    push!(dists, dist)
end

histogram(dists[1], normalize = true, bins = bins, label = L"dim = $2^4$", linewidth = 2, linetype = :stephist)
histogram!(dists[2], normalize = true, bins = bins,label =  L"dim = $2^6$", linewidth = 2, linetype = :stephist)
histogram!(dists[3], normalize = true, bins = bins, label = L"dim = $2^8$", linewidth = 2, linetype = :stephist)
histogram!(dists[4], normalize = true, bins = bins, label = L"dim = $2^{10}$", linewidth = 2, linetype = :stephist)
histogram!(dists[5], normalize = true, bins = bins, label = L"dim = $2^{12}$", linewidth = 2, linetype = :stephist)

#largely dependent on the density

title!("Entanglement spectrum of evolved \n local X vs system size (2 H)")
xlabel!(L"$\lambda$")
savefig("figs/mass_at_zero_finite_size_2_H.png")


#Same thing but now with constant density
ssizes = [2^2, 2^3, 2^4, 2^5]
dists = []

empbins = LinRange(0, 3, 50)

for s in ssizes
    dist = []

    nqubits = convert.(Int,log2(s)*2)

    X = LinearAlgebra.kron(Matrix(I, convert(Integer, s/2), convert(Integer, s/2)), [0 1;1 0], Matrix(I, s, s))

    println((2^6/s)^2)
    for i in 1:convert(Integer, (2^6/s)^2)
        locs = sample(1:nqubits, convert(Int,nqubits/2), replace = false)
        H = insert_op_at_locs(1/sqrt(2)*[1 1 ; 1 -1], nqubits, locs)
        RP1 = Base.permutecols!!(Matrix{Float64}(I, s^2, s^2), Random.randperm(s^2))
        RP2 = Base.permutecols!!(Matrix{Float64}(I, s^2, s^2), Random.randperm(s^2))
        Mat = (RP2*H*RP1)
        #tens = TensorKit.TensorMap(Mat*X*Mat', ℂ^s*ℂ^s, ℂ^s*ℂ^s)
        #U, D, V = tsvd(tens, (1,3), (2,4))
        O = Mat*X*Mat'
        ptranspose!(O, s)
        U,D,V = svd(O)
        map(λ -> push!(dist, λ),D)# diag(D.data))
    end
    push!(dists, dist)
end

histogram(dists[1], normalize = true, bins = empbins, label = L"dim = $2^4$", linewidth = 2, linetype = :stephist)
histogram!(dists[2], normalize = true, bins = empbins,label =  L"dim = $2^6$", linewidth = 2, linetype = :stephist)
histogram!(dists[3], normalize = true, bins = empbins, label = L"dim = $2^8$", linewidth = 2, linetype = :stephist)
histogram!(dists[4], normalize = true, bins = empbins, label = L"dim = $2^{10}$", linewidth = 2, linetype = :stephist)
histogram!(dists[5], normalize = true, bins = bins, label = L"dim = $2^{12}$", linewidth = 2, linetype = :stephist)

pmax = 2^2 #How to explain this?
x = LinRange(0,sqrt(pmax)*.99999, 100)
plot!(x, x->4/(π*pmax)*x*sqrt(pmax/x^2-1), linewidth = 4, label = "Marchenko Pastur", linestyle = :dash, color = :red)


title!("Entanglement spectrum of evolved \n local X vs system size at D_H = 1/2")
xlabel!(L"$\lambda$")
savefig("figs/mass_at_zero_finite_size_fixed_density.png")

#Note: according to "Quantum statistical mechanics of encryption," even permutations are easier to sample from a local
#gate set. Might want to ensure these distributions hold for even permutations!