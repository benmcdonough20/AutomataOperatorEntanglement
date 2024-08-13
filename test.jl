using Pkg; Pkg.activate("automata")
using Plots
using ProgressBars
using StridedViews
using LaTeXStrings
using LinearAlgebra
using StatsBase
using RandomQuantum
using Random

include("tools/pauli_tools.jl")
include("tools/QI_tools.jl")
include("tools/ES_tools.jl")

function read_dat(n, nh)
f = open("/home/ben/Documents/AutonomaProject/data/magical_data/M_N$(n)_nH$(nh)", read = true)
dat = read(f, String)
close(f)
dlist = split(dat, "sample")[2:end]
dlist_split = split.(dlist, "\n")
[parse.(Float64, d[2:end-1]) for d in dlist_split]
end

histogram!(vcat(spacings.(read_dat(10,5))...), normalize = true)

histogram(vcat(read_dat(8,1)...), normalize = true)

histogram(vcat(spacings.(read_dat(10, 5))...), normalize = true)
plot!(LinRange(0,5, 100), x->27/8*(x+x^2)/(1+x+x^2)^(5/2))
xlims!(0,5)

plot!(LinRange(0,2,100), x->1/π*sqrt(4-x^2))




blist = [sum([dlist_nums[30*i + j] for j in 1:30]) for i in 0:1000÷30-1]

bar((1:length(blist))*7/1000*30, blist)
sum(blist)
sum(dlist_nums)

include("tools/imports.jl")

plot

d = 2^6
da = 2^3
db = 2^3
using RandomQuantum
r = ClosedHaarEnsemble(d)
#tesing random overlaps
p = pauli(6, 10)
ret = []
for i in 1:d
    U = rand(r)
    v = U[i,:]
    push!(ret, v'*p*v)
end
std(ret)*d

d = 2^8
da = 2^4 #> db
db = 2^4

dists = []
r = RandomQuantum.GUE(d)
for i in 1:1000
    U = rand(r)
    U .*= db/sqrt(tr(U'U))
    M = sreshape(StridedView(U), (da^2, db^2))
    O = zeros(ComplexF64, da^2, db^2)
    O .= M
    vals = svdvals(O)
    push!(dists, max(vals...)) 
end

dist = [tr(rand(r)) for i in 1:10000]
std(dist)
sqrt(d)/sqrt(2)

dists
histogram(dists, normalize = true)

using Elliptic

β = db^2/da^2
λp = (1+sqrt(β))^2
λm = (1-sqrt(β))^2
xax = LinRange(λm, λp, 200)
xax = LinRange(0.01, 16, 200)
pMP(x) = -(1/π^2)*E((x-16)/x) + (16+x)/(2π^2*x) * K((x-16)/x)
MP(x) = sqrt((x-λm)*(λp-x))/(2π*x*β)

plot!(xax, pMP)

include("tools/pauli_tools.jl")
using LinearAlgebra
using ProgressBars
using StatsBase


d = 2^6
d_a = 2^3
p1 = kron(pauli(3, 10), pauli(3,5))
p2 = kron(pauli(3,10), pauli(3,2))

stds = []
r = RandomQuantum.GUE(d)
for i in 1:10000
    M = rand(r)
    M .*= d_a/sqrt(tr(M'M))
    push!(stds, tr(M))
end
mean(stds)
std(stds)

M[1,1] = 0
heatmap(real.(M))
maximum([abs(M[i,j]) for i in 1:d^2, j in 1:d^2 if i ≠ j])
heatmap(abs.(M) .<= .01)
stds
mean(stds)
var(stds)

n = 10
n2 = 5

sqrt(2^n/tr(M'M))
sqrt(1/2^n)/2


2^(2n2-1) + 2^(n2-1)

d1 = []
d2 = []
for i in ProgressBar(1:128)
    d = 2^(2n2-1) + 2^(n2-1)
    M = randn(d, d)/2^n2 * sqrt(2)
    push!(d1, svdvals(M)...)
    d = 4^(n2) - d
    M = randn(d, d)/2^n2 * sqrt(2)
    push!(d2, svdvals(M)...)
    dist = reverse(sort(dist))
end

histogram(d1, normalize = true)
histogram!(d2, normalize = true)

length(d1)/length(d2)

(2^n2+1)/(2^n2-1)

length(dists)


dist_GOE = []
for i in 1:128
    M = randn(2^n, 2^n)
    M = (M+M')/sqrt(2)
    push!(dist_GOE,  level_spacings(svdvals(M))...)
end

using StatsBase
mean(dist_GOE)
1/2^n2

4^n

length(dist_GOE)

M = randn(2^n, 2^n)
M = (M+M')/sqrt(2)

histogram(dists[dists .<= 25], normalize = true)
histogram(dist_GOE[dist_GOE .<= 25]*8, normalize = true)

d = 2
da = Int(d*(d+1)/2)

level_spacings(dist) = [dist[i]-dist[i+1] for i in 1:length(dist)-1]


lsdist = []
for i in 1:100000
    dist = []
    M1 = randn(d, d)
    M2 = randn(d, d)
    push!(dist, svdvals(M1)...)
    push!(dist, svdvals(M2)...)
    dist = reverse(sort(dist))
    push!(lsdist, level_spacings(dist)...)
end

function rand_PHP(N, nH, nS)
	H = 1/sqrt(2) * [1 1; 1 -1]
    S = [1 0 ; 0 1im]
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	rand_locs = sample(1:N, nS, replace = false)
    Sf = insert_op_at_locs(S, N, rand_locs)
	#Base.permutecols!!(Base.permutecols!!(Hf, randperm(2^N))', randperm(2^N))'
    P1 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P2 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P1 * Hf * Sf * P2
end

histogram!(lsdist[lsdist .<= 25], normalize = true)

dist1 = []
n = 8
X = kron([0 1; 1 0], Matrix(I, 2^(n-1), 2^(n-1)))
for i in 1:100
    O = rand_PHP(8,2, 1)
    M = O*X*O'
    push!(dist1, level_spacings(ES2(M)))
end

#Let's try a couple things to fix this. We can try the trick, we can try inserting S gates, and we can try multiple layers

dist1 = []
n = 8
X = [0 1; 1 0]
#P = kron(exp(1im * 5π/16 * X)*[1 0; 0 -1]*exp(-1im * 5π/16 * X), Matrix(I, 2^(n-1), 2^(n-1)))
P = kron([1 0; 0 -1], Matrix(I, 2^(n-1), 2^(n-1)))

for i in 1:100
    O = rand_PHP(8,3,3)#rand_PHP(8,1,1) * rand_PHP(8,1,1) * rand_PHP(8,1,1)
    M = O*P*O'
    push!(dist1, ES2(M))
end

histogram!(vcat(dist1...), normalize = true)
histogram!(vcat(dist2...), normalize = true)

plot!(LinRange(0,2,100), x->1/π*sqrt(4-x^2))
xlims!(0,5)
xax = LinRange(0,2, 100)
Zb = 8/27
b = 1
plot!(xax, r-> 1/Zb * (r+r^2)^b/(1 + r + r^2)^(3/2*b + 1), color = :black, linestyle = :dash, linewidth = 3, label = "WD, β = 1")

f = open("M_N12_nH1.dat", "r")
r = read(f, String)
close(f)

dat = split(r, ",")
dists = [parse.(Float64,split(d, "\n")[2:end-1]) for d in dat[2:end]]

histogram(vcat(dists...))

using Plots
using Arrow

df = Arrow.Table("X_N6_nH1.arrow")

histogram(vcat(df...), normalize = true)
plot!(LinRange(0, 2, 100), x->1/π*sqrt(4-x^2))

lambda = []
for i in ProgressBar(1:10)
    O = rand_PHP(10, 0, 1)
    X = kron([0 1; 1 0], Matrix(I, 2^9, 2^9))

    M = O*X*O'


    push!(lambda, ES2(M)...)
end

histogram(lambda)
cm = countmap(round.(lambda, digits = 14))

bar(collect(keys(cm)), collect(values(cm)))

ns = [1,2,3,4]
dists = []
for n in ns
    r = RandomQuantum.GUE(2^(2*n))
    dist = []
    for i in 1:500
        M = rand(r)
        ptranspose!(M)
        vals = svdvals(M)
        vals = vals.^2 ./ norm(vals)^2
        push!(dist, -sum(vals .* log.(vals)))
    end
    push!(dists, mean(dist))
end

plot(ns, dists)
plot!(LinRange(0, 4, 10), x->x*log(4)-1/2)

function rand_PHP(N, nH)
	H = 1/sqrt(2) * [1 1; 1 -1]
	#Rx = exp(1im * π/4 * [0 1; 1 0])
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	rand_loc = rand(1:N)
    #Rxf = insert_op_at_locs(Rx, N, [rand_loc])
    P1 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P2 = Base.permutecols!!(Matrix(I, 2^N, 2^N), randperm(2^N))
    P1 * Hf * P2
end

O = rand_PHP(10, 1)
O == O'

X = kron([0 1; 1 0], Matrix(I, 2^7, 2^7))
dist = []
for i in 1:100
    Mat = rand_PHP(8, 2)
    O = Mat * X * Mat'
    vals = ES2(O)
    push!(dist, vals...)
end

histogram!(dist, normalize = true)


O = rand_PHP(10, 10)
X = kron([0 1; 1 0], Matrix(I, 2^9, 2^9))
M = O*X*O

histogram([M[i,j] for i in 1:2^10 for j in 1:2^10])
xlims!(-.1, 1)

dists = []
for s in [2, 3, 4, 5]
    dist = []
    for i in 1:100
        M = Matrix{Float64}(abs.(randn(2^(2s), 2^(2s))))
        vals = ES2(M)
        vals = vals.^2 / norm(vals)^2
        push!(dist, S(vals))
    end
    push!(dists, mean(dist))
end

dists
s = 4
M = Matrix{Float64}(abs.(randn(2^(2s), 2^(2s))))
vals = ES2(M)
vals = vals.^2 / norm(vals)^2
histogram(vals)
xlims!(0, .1)
plot(dists)
plot!(2:5, x->(1-μ2/α)*x*log(4))

μ2 = mean(abs.(randn(2^(2s), 2^(2s)))).^2
α = mean(abs.(randn(2^(2s), 2^(2s))).^2)


#string wavefunction demonstration

n = 6
plots = []
nh = 3
    O = rand_PHP(6, nh)
    X = kron([0 1; 1 0], Matrix(I, 2^(n-1), 2^(n-1)))

    realidxs = [i for i in 0:2^n-1 if isreal(pauli(Int(n/2), i))]
    cmplidxs = [i for i in 0:2^n-1 if !isreal(pauli(Int(n/2), i))]

    Mr = [tr(kron(pauli(Int(n/2), realidxs[i]), pauli(Int(n/2),realidxs[j]))*(O*X*O')) for i in 1:length(realidxs), j in 1:length(realidxs)]
    Mc = [tr(kron(pauli(Int(n/2), cmplidxs[i]), pauli(Int(n/2),cmplidxs[j]))*(O*X*O')) for i in 1:length(cmplidxs), j in 1:length(cmplidxs)]

    p1 = heatmap(real.(Mr))
    p2 = heatmap(real.(Mc))

    maximum(real.(Mr))
    maximum(real.(Mc))
    p = plot(p1, p2, layout = (2,1), size = (300,500))
    push!(plots, p)

plots[1]
plots[2]
plots[3]

plot(plots..., layout = (1, 3), size = (1100, 600), plot_title = L"String wavefunction $n_H = 0$ $n_H = 1$ $n_H = 6$ $\Re$ $\Im$")
display(p)
savefig("figs/Pauli_wf.svg")

n = 10
nh = 2
    f = open("/home/ben/Documents/AutonomaProject/data/magical_data/Z_N$(n)_nH$(nh)", "r")
r = read(f, String)

dat = split(r, "sample")
datparse = [parse.(Float64, split(d, "\n")[2:end-1]) for d in dat]

histogram(vcat(datparse...), normalize = true)

dists = []
Z = kron([1 0; 0 -1], Matrix(I, 2^7, 2^7))
for i in 1:10*4
    Mat = rand_PHP(8,4)
    O = Mat*Z*Mat'

    vals = ES2(O)
    push!(dists, vals)
end
histogram!(vcat(dists...), normalize = true)
plot!(LinRange(0,2,100), x->sqrt(4-x^2)*1/π)
d = vcat(spacings.(dists)...)
histogram(d[d.<=25])
xlims!(0, 5)

spacs1 = []
d = RandomQuantum.GUE(2)
for i in 1:100*2^12
    re = randn()
    m,_ = eigen(rand(d))
    s = sort([m; [re]])
    push!(spacs1, (s[3]-s[2])/(s[2]-s[1]))
end

histogram(spacs1[spacs1 .< 100], normalize = true, linetype = :stephist)
xlims!(0, 10)
plot!(LinRange(0, 10, 100), x-> 3/4*(1+x)/(1+x+x^2)^(3/2), linewidth = 3, color = :red)

spacs = []
for i in 1:1000
    M = randn(2^6, 2^6)
    M = M+M'
    push!(spacs, spacings(ES2(M))...)
end
histogram!(spacs[spacs .< 100], normalize = true, linetype = :stephist)
a = 3*sqrt(3)/(2*π)
b = 1
plot!(LinRange(0, 10, 100), x-> a/(1+(1+x/b)^2 - 1/3*(2+x/b)^2), linewidth = 3, color = :red)
xlims!(0,4)


d = 2000
P = rand(d,d) .< 1/d
M = P*P'

1/d*tr(M)
1/d*tr(M^2)
1/d*tr(M^3)
1/d*tr(M^4)
1/d*tr(M^5)

k = 2
iter1 = Base.product([1:d for i in 1:k]...)
iter2 = Base.product([1:d for i in 1:k]...)

f = open("MPOs/X_data_14.dat", "r")

d = read(f, String)
arrs_str = split(d, "\n")
arrs_str = [a[2:end-1] for a in arrs_str]
arrs = split.(arrs_str, ",")
arrs = [parse.(Float64, a) for a in arrs[1:end-1]]

using Plots
histogram(arrs[1])
histogram!(arrs[2])
histogram!(arrs[3])

sum(arrs[1].^2)/2^14
sum(arrs[2].^2)/2^14

mean(arrs[1].^4)
mean(arrs[2].^6)

mean([mean(arrs[i].^10) for i in 1:256])
mean([mean(arrs[i].^8) for i in 1:256])
mean([mean(arrs[i].^6) for i in 1:256])
mean([mean(arrs[i].^4) for i in 1:256])

dist = []
r = RandomQuantum.ClosedHaarEnsemble(2^8)
X = kron(Matrix(I, 2^3, 2^3), [0 1; 1 0], Matrix(I, 2^4, 2^4))
for i in 1:200
    M = rand(r)
    O = M'*X*M
    ptranspose!(O)
    push!(dist, svdvals(O))
end
histogram(vcat(dist...))

spectralformfactor(dist, t) = abs(sum(exp.(1im .* dist * t)))^2

plot(10 .^ LinRange(-2, 10, 1000), t->mean(spectralformfactor.(dist, t))*2^6, xaxis = :log, yaxis = :log, label = L"U^\dagger X U")

plot(10 .^ LinRange(-2, 10, 1000), t->mean(spectralformfactor.(arrs, t)), xaxis = :log, yaxis = :log, label = L"P^t X P")
plot!(10 .^ LinRange(-2, 10, 1000), t->spectralformfactor(arrs[1], t), xaxis = :log, yaxis = :log, label = L"P^t X P")

xlabel!(L"\theta")
ylabel!(L"\langle g_M(\tau) \rangle")
title!("Spectral form factor")
savefig("figs/spectral_form_factor.png")

dist = randn(10000)

plot(10 .^ LinRange(-2, 10, 1000), t->mean(spectralformfactor(dist, t)), xaxis = :log, yaxis = :log, label = L"U^\dagger X U")
