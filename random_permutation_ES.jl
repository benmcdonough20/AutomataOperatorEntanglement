include("tools/imports.jl")
include("plotting_defaults.jl")
include("tools/ES_tools.jl")

using Kronecker
using LaTeXStrings

#Workflow for getting a random permutation, writing its cycle decomposition, and its entanglement spectrum
Da = 6
Db = 6

#average permutation spectrum Bernoulli vs Permutation
samples = 2^5
N = 2^(Da+Db)
ss = [6,8,10]
dists = []

dist_bern = []
for i in ProgressBar(1:samples)
    P = convert.(Float64,rand(N,N) .<= 1/N)
    ptranspose!(P)
    push!(dist_bern, svdvals(P)...)
end

dist_perm = []
for i in ProgressBar(1:samples)
    P = rand_perm_mat(N);
    ptranspose!(P)
    push!(dist_perm, svdvals(P)...)
end

mp1 = countmap(round.(dist_bern, digits = 13))
mp2 = countmap(round.(dist_perm, digits = 13))
d = 0
ks = union(keys(mp1), keys(mp2))
for k in ks
    x = haskey(mp1, k) ? mp1[k] : 0
    y = haskey(mp2, k) ? mp2[k] : 0
    d+= abs(x-y)
end
push!(dists, d/length(ks))


# plotting a discrete distribution
mp1 = countmap(round.(dist_bern, digits = 13))
mp2 = countmap(round.(dist_perm, digits = 13))
bar(collect(keys(mp1)), collect(values(mp1))./sum(values(mp1)), yaxis = :log, linecolor = :blue, alpha = .8, color = :blue, label = "permutations")
bar!(collect(keys(mp2)), collect(values(mp2))./sum(values(mp2)), yaxis = :log, linecolor = :green, alpha = .8, color = :green, label = "dilute graphs")
est = exp(-exp(-1)) + 2*exp(-1)-1
hline!([est], label = L"e^{-1/e}+2/e-1", color = :red)
ylims!(.5/length(dist_bern), 1)
title!("Random permutation ES vs spectrum of dilute random graph")
xlabel!(L"\lambda")
ylabel!(L"P(\{\lambda\})")
savefig("figs/paper/bern_vs_perm.svg")

#average permutation spectrum Bernoulli vs Permutation
samples = 2^18
ss = [6,8,10, 12]
dists = []
for s in ss
    N = 2^s
    dist_bern = []
    for i in ProgressBar(1:samples*2.0^(-s))
        P = convert.(Float64,rand(N,N) .<= 1/N)
        ptranspose!(P)
        push!(dist_bern, svdvals(P)...)
    end

    dist_perm = []
    for i in ProgressBar(1:samples * 2.0^(-s))
        P = rand_perm_mat(N);
        ptranspose!(P)
        push!(dist_perm, svdvals(P)...)
    end

    mp1 = countmap(round.(dist_bern, digits = 13))
    mp2 = countmap(round.(dist_perm, digits = 13))
    d = []
    ks = union(keys(mp1), keys(mp2))
    for k in ks
        x = haskey(mp1, k) ? mp1[k]/length(dist_bern) : 0
        y = haskey(mp2, k) ? mp2[k]/length(dist_bern) : 0
        push!(d, abs(x-y))
    end
    push!(dists, maximum(d))
end

dists
plot(6:2:12, log10.(dists), label = "", color = :black, marker = :o)
xticks!(6:2:12)
xlabel!(L"N")
ylabel!(L"\ln \Vert ES(P_N)-ES(B_N) \Vert_{\infty}")
title!("Convergence of automaton ES and dilute graph spectral measure")
savefig("figs/paper/perm_bern_convergence_linplot.svg")

mp2[0]/length(dist_perm)

#Workflow for getting a random permutation, writing its cycle decomposition, and its entanglement spectrum
Da = 6
Db = 6

#average permutation spectrum Bernoulli vs Permutation
samples = 10
N = 2^(Da+Db)
dist_one = []
samples = 1
P = convert.(Float64,rand(N,N) .<= 1/N)
ptranspose!(P)
dist_one = svdvals(P)

Da = 4
Db = 4
dist_many = []
samples = 200
for i in ProgressBar(1:samples)
    P = rand_perm_mat(2^Da*2^Db);
    ptranspose!(P)
    push!(dist_many, svdvals(P)...)
end

# plotting a discrete distribution
mp1 = countmap(round.(dist_one, digits = 13))
mp2 = countmap(round.(dist_many, digits = 13))
bar(collect(keys(mp1)), collect(values(mp1))./sum(values(mp1)), yaxis = :log, linecolor = :blue, alpha = .8, color = :blue, label = L"N = 12 (\times 1)")
bar!(collect(keys(mp2)), collect(values(mp2))./sum(values(mp2)), yaxis = :log, linecolor = :green, alpha = .8, color = :green, label = L"N = 8 (\times 256)")

hline!([est], label = L"e^{-1/e}+2/e-1", color = :red)
ylims!(.5/length(dist_one), 1)
xlabel!(L"\lambda")
ylabel!(L"P(\{\lambda\})")
title!("Concentration of automaton ES")
savefig("figs/paper/perm_ES_concentration.svg")


empbins = LinRange(0, 3, 266)
histogram(dist_one, normalize = true, bins = empbins, alpha = .5)
histogram!(dist_many, normalize = true, bins = empbins, alpha = .5)

mp1[0]/length(dist_one)
mp2[0]/length(dist_many)

est = exp(-exp(-1)) + 2*exp(-1)-1
hline!([est], label = L"e^{-1/e}+2/e-1", color = :red)
ylims!(.5/length(dist_bern), 1)
title!("Random permutation ES vs spectrum of sparse random graph")
xlabel!(L"\lambda")
ylabel!(L"P(\{\lambda\})")
savefig("figs/paper/bern_vs_perm.svg")

mp2[0]/length(dist_perm)

Da = 4
Db = 4
D = Da+Db
Z = kron([1 0; 0 -1], Matrix(I, 2^(D-1), 2^(D-1)))
X = kron([0 1; 1 0], Matrix(I, 2^(D-1), 2^(D-1)))
Y = kron([0 1im; -1im 0], Matrix(I, 2^(D-1), 2^(D-1)))
op = kron([0 0; 1 0], Matrix(I, 2^(D-1), 2^(D-1)))
dist_Z = []
samples = 10
for i in ProgressBar(1:samples)
    P = rand_perm_mat(2^Da*2^Db);
    M = 1/sqrt(2^Db)*P*(Z)*P' #same for Y and X -- for all off-diagonal operators?
    ptranspose!(M)
    push!(dist_Z, svdvals(M)...)
end


op = kron([0 0; 1 0], Matrix(I, 2^(D-1), 2^(D-1)))
dist_1 = []
dist_2 = []
samples = 1000
for i in ProgressBar(1:samples)
    P = rand_perm_mat(2^Da*2^Db);
    M = P*(X)*P' #same for Y and X -- for all off-diagonal operators?
    ptranspose!(M)
    push!(dist_1, svdvals(M)...)
    M = P*(Z+X)*P' #same for Y and X -- for all off-diagonal operators?
    ptranspose!(M)
    push!(dist_2, svdvals(M)...)
end

dist_1_rounded = round.(dist_1, digits = 13)
dist_2_rounded = round.(dist_2, digits = 13)

cm1 = countmap(dist_1_rounded)
cm2 = countmap(dist_2_rounded)
cm1[0]/length(dist_1)
cm2[0]/length(dist_1)
bar(collect(keys(cm1)), collect(values(cm1)))
bar!(collect(keys(cm2)), collect(values(cm2)))

histogram(dist_1_rounded)
histogram!(dist_2_rounded)

#Why not just restrict ourselves to Hermitian obserables??

Da = 2
Db = 2
D = 4

P = rand_perm_mat(2^Da*2^Db);
M = P*(op+op')*P' #same for Y and X -- for all off-diagonal operators?

M = P*(op)*P' #same for Y and X -- for all off-diagonal operators?

dist_rounded = round.(dist_Z, digits = 5)
cm = countmap(dist_rounded[dist_rounded .!= 0])
bar(collect(keys(cm)), collect(values(cm)))
xlims!(-.1,2)


P = rand_perm_mat(2^Da*2^Db);
M = P*(Z)*P' #same for Y and X -- for all off-diagonal operators?

M = diagm(rand([-1, 1], 2^D))

ptranspose!(M)
vals = svdvals(M)

dist_rounded = round.(vals, digits = 5)
cm = countmap(dist_rounded[])
bar!(collect(keys(cm)), collect(values(cm)))

Da = 2
Db = 2
D = Da+Db

Z = kron([1 0; 0 -1], Matrix(I, 2^(D-1), 2^(D-1)))
X = kron([0 1; 1 0], Matrix(I, 2^(D-1), 2^(D-1)))
Y = kron([0 1; -1m 0], Matrix(I, 2^(D-1), 2^(D-1)))

P = rand_perm_mat(2^Da*2^Db);

P*Z*P'
P*X*P
M = P*Y*P'
ptranspose!(M)
M*M'

Da = 6
Db = 6
D = 12
dist = []
for i in ProgressBar(1:100)
ψ = zeros(Float64, 2^D)
ψ[rand(1:2^D)] = 1
P = rand_perm_mat(2^Da*2^Db);
H = kron([1/sqrt(2)*[1 1; 1 -1] for i in 1:D]...)
ψ = H*ψ
ψ = P*ψ
ψ = H*ψ
v = zeros(ComplexF64, 2^Da, 2^Db)
v .= sreshape(StridedView(ψ), (2^Da, 2^Db))
v .*= sqrt(2^Da)
vals = svdvals(v)
push!(dist, vals...)
end

histogram(dist, normalize = true)
plot!(LinRange(0, 2, 100), x->1/π * sqrt(4-x^2))


samples = 2^16
ss = [6,8,10,12]
atoms = []
stds = []
for s in ss
    N = 2^s
    dist_perm = []
    for i in ProgressBar(1:samples * 2.0^(-s))
        #P = rand_perm_mat(N);
        #Z,M = gen_large_cluster(2^s)
        M = rand(2^s, 2^s) .<= 1/2^s
        push!(dist_perm, svdvals(M))
    end
    ests = []
    for i in 1:length(dist_perm)
        samp = dist_perm[i]
        mp = countmap(round.(samp, digits = 10))
        push!(ests, sum([mp[k] for k in collect(keys(mp)) if mp[k]>2])/length(samp))
    end
    mp = countmap(round.(vcat(dist_perm...), digits = 10))
    push!(atoms, sum([mp[k] for k in collect(keys(mp)) if mp[k]>2])/length(vcat(dist_perm...)))
    push!(stds, std(ests)/sqrt(length(dist_perm)))
end

plot(6:2:12, atoms, linewidth = 2, yerr = stds)
title!("Size of singular component of spectrum")
xlabel!("N")
ylabel!("p")

plot!(LinRange(6, 12, 100), x->1/sqrt((x-6)))
savefig("figs/paper/singularfrac.svg")

moms = Dict()
sizes = [6,8,10,12, 14]
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

f = open("rand_perm_moms.dat", "w")
println(f, moms)
close(f)

for s in sizes
    for i in 1:4
        moms[s][i] = real(moms[s][i])
    end
end

f = open("MPOs/X_data_14.dat", "r")

d = read(f, String)
arrs_str = split(d, "\n")
arrs_str = [a[2:end-1] for a in arrs_str]
arrs = split.(arrs_str, ",")
arrs = [parse.(Float64, a) for a in arrs[1:end-1]]
close(f)

moms[14] = [[],[],[], []]
for a in arrs
    for i in 1:4
        push!(moms[14][i], sum(a .^ ((i+1)*2)) / 2^14)
    end
end

include("tools/plotting_defaults.jl")

p = plot(guidefontsize = 14)
moms_ideal = [3, 12, 57, 303]
for (k,m) in zip(1:4, moms_ideal)
    plot!(sizes, [(m-mean(moms[s][k]))/m for s in sizes], yerr = [std(moms[s][k])/(m*sqrt(length(moms[s][k]))) for s in sizes], yaxis = :log, label = "k = $(k+1)", color = seaborn[k], linewidth = 2, marker = "x")
end
display(p)

xlabel!(L"N")
ylabel!(L"\frac{\langle \lambda^k \rangle - m_k}{m_k}")
title!(L"Convergence of $X(\tau)$ OES moments to Bernoulli ensemble")

savefig("figs/paper/moments_convergence.svg")

log10((303-mean(moms[14][4]))/303)