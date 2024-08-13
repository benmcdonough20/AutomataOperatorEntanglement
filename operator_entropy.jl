include("tools/imports.jl")
import Base.rand

#check that this is approximately 1--proper normalization

sizes = [2^4, 2^6, 2^8, 2^10]
samples = 100
avgs_Haar = []
avgs_GUE = []
for s in ProgressBar(sizes)
    avg1 = 0
    avg2 = 0
    dim = round(Int, sqrt(s))
    d = RandomQuantum.ClosedHaarEnsemble(s)
    q = RandomQuantum.GUE(s)
    for i in 1:samples
        M = rand(d)
        M ./= sqrt(tr(M'M))
        P = rand(q)
        P ./= sqrt(tr(P'P))
        avg1 += S(ES(M, dim))
        avg2 += S(ES(P, dim))
    end
    push!(avgs_Haar, avg1/samples)
    push!(avgs_GUE, avg2/samples)
end


plot(0:1:5, x->x*log(2)-1/2, label = L"L\ln(2)-1/2", color = :black, linewidth =2, linestyle = :dash)
plot!(0:1:5, x->2*x*log(2)-1/2, label = L"2L\ln(2)-1/2", color = :black, linewidth = 2, linestyle = :dash)
plot!(0:1:5, x->2*x*log(2)-0.86, label = L"2L\ln(2)-ΔS_{P}", color = :black, linewidth = 1)
scatter!(2:1:5, avgs_Haar, label = "Haar", legend = :topleft)
scatter!(2:1:5, avgs_GUE, marker = :xcross, markersize = 5, label = "GUE")
xlabel!(L"$L$")
ylabel!(L"$S_{op}$")
scatter!(2:1:5, avgsX, label = L"$\sigma X\sigma^\dagger$")
scatter!(2:1:5, avgsZ, label = L"$\sigma Z \sigma^\dagger$")

avgsX[1] - 2*2*log(2)

title!("Operator entanglement entropy vs system size")
savefig("figs/paper/operator_entanglement.svg")

struct rand_conj_perm_dist
    op
    dim
end

import Base.rand

function rand(dist::rand_conj_perm_dist)
    D = dist.dim
    P = Base.permutecols!!(Matrix{Float64}(I, D, D), Random.randperm(D))
    P*dist.op*P'
end

#Does the entanglement entropy saturate for a permutation in the same way?
sizes = [2^4, 2^6, 2^8, 2^10]
samples = 100

X = [0 1; 1 0]
Z = [1 0; 0 -1]
#ops = [sqrt(2)*(a*X + (1-a)*Z)/norm(a*X+(1-a)*Z) for a in LinRange(0, 1, 10)]
ops = [X, Z]
avgss = [[] for op in ops]
avgs_GUE = []
avgs_Haar = []

for s in sizes
    dim = round(Int, sqrt(s))
    for (k,loc_op) in enumerate(ops)
        avg = 0
        op = kron(
            Matrix{Float64}(I, dim÷2,dim÷2),
            loc_op,
            Matrix{Float64}(I, dim, dim)
        )
        d = rand_conj_perm_dist(op, s)
        for i in ProgressBar(1:samples)
            M = rand(d)/sqrt(s)
            avg += S(ES(M, dim))
        end
        push!(avgss[k], avg/samples)
    end
end

for s in sizes
    dHaar = RandomQuantum.ClosedHaarEnsemble(s)
    dGUE = RandomQuantum.GUE(s)
    avg_Haar = 0
    avg_GUE = 0
    dim = round(Int, sqrt(s))
    for i in ProgressBar(1:samples)
        MHaar = rand(dHaar)/sqrt(s)
        MGUE = sqrt(2)*rand(dGUE)/s
        avg_Haar += S(ES(MHaar, dim, dim)) 
        avg_GUE += S(ES(MGUE, dim, dim)) 
    end
    push!(avgs_Haar, avg_Haar/samples)
    push!(avgs_GUE, avg_GUE/samples)
end

avgs_perm = []
for s in sizes
    avg_perm = 0
    dim = round(Int, sqrt(s))
    for i in ProgressBar(1:samples)
        Mperm = rand_perm_mat(s)/sqrt(s)
        avg_perm += S(ES(Mperm,dim, dim)) 
    end
    push!(avgs_perm, avg_perm/samples)
end


cols = [:red, :green, :blue, :orange, :cyan ,:purple]
cy = Base.Iterators.cycle(cols)

slopes = []
p = plot()
x = [i for i in 1:5]
for (c, k,avgs) in zip(cy,LinRange(0, 1, 10), avgss)
    scatter!(avgs, color = c, label = "")
    #a = x[3:end] ⋅ (avgs[3:end].+1/2) / sum(x[3:end].^2)
    b,a = linear_fit(x[2:end], avgs[2:end])
    push!(slopes, a)
    #plot!(0:1:5, x->a*x-1/2, label = "", color = c)
    #plot!(0:1:5, x->(sqrt(2)*k/norm(X*k+(1-k)*Z)+1)log(2)*x-1/2, label = "", color = c)
end

scatter(avgss[1], color = :blue, label = "X", dpi = 300)
scatter!(avgss[2], color = :red, label = "Z")
scatter!(avgs_Haar, color = :green, label = "Haar/GUE")
#scatter!(avgs_perm, color = :purple, label = "perm")

plot!(0:1:5, x->2*x*log(2)-1/2, label = "Random operator entanglement", color = :black, linewidth = 3, linestyle = :dash, dpi = 300)
plot!(0:1:5, x->2*x*log(2)-1, label = "Random operator entanglement", color = :black, linewidth = 3, linestyle = :dash, dpi = 300)
plot!(0:1:5, x->x*log(2)-1/2, label = "Random state entanglement", color = :black, linewidth = 3, linestyle = :dash)

display(p)
xlabel!(L"$L$")
ylabel!(L"$S_{op}$")
#title!(L"Interpolation of $S_{op}$ between $X$ and $Z$")
title!(L"Dependency of $S_{op}$ on initial operator")
savefig("figs/operator_entanglement_vs_system_size.png")

plot(LinRange(0,1,10),[(s-log(2))*norm(k*X+(1-k)*Z) for (s,k) in zip(slopes, LinRange(0,1,10))], color = :red, linewidth = 2, label = "fit")
plot!(LinRange(0,1,10),LinRange(0,1,10), color = :green, label = "linear interpolation")
xlabel!("α")
ylabel!(L"$(y-1/2)/ln(2)-1$")
title!("Operator entanglement entropy slope interpolation \n between X and Z")
savefig("figs/slope_interpolation.png")
slopes

rand_perm_mat = N->Base.permutecols!!(Matrix{Float64}(I, N, N), randperm(N))

eigs = []
for i in ProgressBar(1:50)
    Perm = rand_perm_mat(2^10)
    map(x->push!(eigs, x),ES(Perm, 2^5, 2^5))
end

bins = LinRange(0,6,72)
histogram(eigs, normalize = true, bins = bins, linetype = :stephist, label = "", linewidth = 2, color = :blue, dpi = 300)
ylabel!("PDF")
xlabel!(L"$\lambda^2$")
title!("Random permutation ES")
savefig("figs/rand_perm_ES")

samples = 100
avgsX = []
avgsZ = []

for s in sizes
    avgX = 0
    avgZ = 0
    dim = round(Int, sqrt(s))
    X = kron(
        Matrix{Float64}(I, dim÷2,dim÷2),
        [0 -1im; 1im 0],
        Matrix{Float64}(I, dim, dim)
    )
    Z = kron(
        Matrix{Float64}(I, dim÷2,dim÷2),
        [1 0; 0 -1],
        Matrix{Float64}(I, dim, dim)
    )
    dX = rand_conj_perm_dist(X, s)
    dZ = rand_conj_perm_dist(Z, s)
    for i in ProgressBar(1:samples)
        MX = rand(dX)/sqrt(s)
        MZ = rand(dZ)/sqrt(s)
        #M = Base.permutecols!!(Matrix{Float64}(I, s, s), Random.randperm(s))/sqrt(s)
        avgX += S(ES(MX, dim))
        avgZ += S(ES(MZ, dim))
    end
    push!(avgsX, avgX/samples)
    push!(avgsZ, avgZ/samples)
end