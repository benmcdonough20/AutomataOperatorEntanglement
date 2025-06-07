include("../tools/imports.jl")

#Phases plot
rel_phases = π .+(0:8)*1/16*π
dists = []
n=10

d = 2^n
r = GOE(d)

for ϕ in ProgressBar(rel_phases)
    U = exp(1im*ϕ*[0 1; 1 0])
    X = U*[1 0; 0 -1]*U'
    X = kron(X, Matrix(I, 2^(n-1), 2^(n-1)))
    dist = []
    for i in ProgressBar(1:100)
        O = rand(r)
        M = O*X*O'
        vals = ES2(M)
        push!(dist, spacings(vals)...)
    end
    push!(dists, dist)
end

p = plot()
labels = ["$(i)π/16" for i in 16:(16+8)]
labels[1] = "π"
for (i,dist) in enumerate(dists)
    if (i-1)%4 == 0
        c = seaborn[1]
    else
        c = seaborn[4]
    end
    lab = ""
    if i == 1
        lab = L"ϕ = n\pi/4"
    end
    if i == 2
        lab = L"ϕ \neq n\pi/4"
    end
    h = bin(dist, 0, 4, 50)
    xax = LinRange(0, 4, 50) .+ 4/100
    scatter!(xax, h, color = c, label = lab, markerstrokewidth = 0, marker = :o, markersize = 3)
end

xlims!(0, 4)
xlabel!(L"r")
ylabel!(L"p(r)")
title!("Rotation angle vs ES spacing ratios")

xax = LinRange(0,4, 100)
Zb = 8/27
b = 1
plot!(xax, r-> 1/Zb * (r+r^2)^b/(1 + r + r^2)^(3/2*b + 1), color = seaborn[4], label = "WD, β = 1")
plot!(xax, x-> 3/4*(x+1)/(1+x+x^2)^(3/2), color = seaborn[1], label = "Surmise for GOE")

savefig("figures/fig11.svg")