#Fig 3: Entanglement spectrum of unitary-evolved operator at different cut sizes
include("../tools/imports.jl")
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