include("../tools/imports.jl")

# Fig. 2: Entanglement level spacing ratios in random observables and comparison to classical ensembles
δWD(x, C, cb, β) = C/(1+x)^2*((x+1/x)^(-β)-cb*(x+1/x)^(-(β+1)))

function WD(β)
    if β == 0 #For real matrices
        return x->3/4*(x+1)/(1+x+x^2)^(3/2)
    elseif β == 1
        C = 0.233378
        cb = 2*(π-2)/(4-π)
        return x->27/8*(x+x^2)/(1+x+x^2)^(5/2)+δWD(x, C, cb, β)
    elseif β == 2
        C = 0.578846
        cb = 4*(4-π)/(3π-8)
        return x->81/4*sqrt(3)/π*(x+x^2)^2/(1+x+x^2)^(4) + δWD(x, C, cb, β)
    end
end

## Make data
d = 2^10
n = Int(log2(d))
samples = 250

rngs = [()->randUnitary(d), ()->rand(GUE(d)), ()->rand(GOE(d))]

dists_classical = [[] for r in rngs]
for i in ProgressBar(1:samples)
    for (i,r) in enumerate(rngs)
        M = convert(Matrix{ComplexF64}, r())
        M .*= sqrt(d/tr(M*M'))
        push!(dists_classical[i], ES2(M))
    end
end

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

## Make plot
xend = 3
xax = LinRange(0, xend, 100)

plot(xax, WD(2), color = seaborn[1], label = L"WD, $\beta = 2$")
plot!(xax, WD(1), color = seaborn[2], label = L"WD, $\beta = 1$")
plot!(xax, WD(0), label = L"\text{Surmise for GOE}", color = seaborn[3])

labels = ["Haar", "GUE", "GOE"]

for (i,dist) in enumerate(dists_classical)
    nbins = 60
    x = LinRange(0,xend,nbins)[1:end-1]
    spacdist = vcat(spacings.(dist)...)
    histplot = [sum(x[i-1] .<= spacdist .< x[i])*(nbins/xend)/length(spacdist) for i in 2:length(x)]
    scatter!(x[1:end-1] .+ xend/(2*nbins), histplot, markershape = :+, color = seaborn[i], label = labels[i], markersize = 7)
end

for (i,dist) in enumerate(dists_evo)
    nbins = 60
    x = LinRange(0,xend,nbins)[1:end-1]
    spacdist = vcat(spacings.(dist)...)
    histplot = [sum(x[i-1] .<= spacdist .< x[i])*(nbins/xend)/length(spacdist) for i in 2:length(x)]
    scatter!(x[1:end-1] .+ xend/(2*nbins), histplot, marker = :o, color = seaborn[i+1], label = labels[i], markerstrokewidth =0, markersize = 4)
end

xlabel!(L"r")
ylabel!(L"p(r)")
savefig("figs/paper/ratios_comparison.svg")

#inset
xend = 2.2
plot(guidefontsize = 16, axeslabelfontsize = 16, tickfontsize = 16)
for (i,dist) in enumerate(dists_classical)
    nbins = 60
    x = LinRange(0,xend,nbins)[1:end-1]
    plot!(x .+ xend/(2*nbins), bin(vcat(dist...), 0, xend, nbins), palette = seaborn, linewidth = 5, label = "")
end

for (i,dist) in enumerate(dists_evo)
    nbins = 120
    x = LinRange(0,xend,nbins)[1:end-1]
    plot!(x .+ xend/(2*nbins), bin(vcat(dist...), 0, xend, nbins), palette = seaborn, linewidth = 5, label = "")
end

plot!(LinRange(0, 2, 100), x->1/π*sqrt(4-x^2), color = :black, linestyle=:dash, linewidth = 3, label = "semicircle")

xlabel!(L"\sqrt{\lambda}")
ylabel!(L"p(\sqrt{\lambda})")

#Fig 3: Entanglement spectrum of unitary-evolved operator at different cut sizes
dist_evos = []
samples = 1
for i in ProgressBar(5:7)
    dist_evo = []
    Da = 2^i
    Db = 2^(10-i)
    X = LinearAlgebra.kron(Matrix(I, convert(Integer, Da/2), convert(Integer, Da/2)), [0.0 1.0;1.0 0.0], Matrix(I, Db, Db))
    X .*= Db / sqrt(tr(X'*X)) #normalization
    d = RandomQuantum.ClosedHaarEnsemble(Da*Db)
    for i in ProgressBar(1:samples)
        RU = rand(d)
        M = RU*X*RU'
        tens = TensorKit.TensorMap(M, ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db)
        U, D, V = tsvd(tens, (1,3),(2,4))
        push!(dist_evo, (D.data.^2)...)
    end
    push!(dist_evos, dist_evo)
end

plot()
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

    dist_binned = [sum(emp_bins[i] .<= dist .< emp_bins[i+1]) for i in 1:length(emp_bins)-1] / (length(dist)*δ)
    scatter!(emp_bins[1:end-1] .+ δ/2, dist_binned, m = :x, markersize = 5, color = seaborn[i%3+1], label = "n_A = $(i)")
    MP(x) = 1/(2π*x*β) * sqrt((λp-x)*(x-λm))
    plot!(xax, MP, linestyle = :dash, linewidth = 2, color = seaborn[i%3+1], label = "")
end

title!(L"Concentration of ES of $UXU^\dagger$ to MP")
xlabel!(L"\lambda")
ylabel!(L"p(\lambda)")