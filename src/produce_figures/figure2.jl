# Fig. 2: Entanglement level spacing ratios in random observables and comparison to classical ensembles
include("../tools/imports.jl")

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
savefig("../../figures/figure2.svg")

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
savefig("../../figures/figure2_inset.svg")