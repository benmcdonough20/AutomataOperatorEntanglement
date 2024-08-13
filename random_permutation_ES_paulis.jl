#another approach: examine the distribution of vector coefficients in Pauli basis
n = 5
D = 2^n
dist_XZ = []
dist_X = []
for i in ProgressBar(1:1000)
    M = rand_perm(D)
    for j in 0:D^2-1
        if j < D
            push!(dist_X, tr(pauli(n, j)'*M))
        else
            push!(dist_XZ, tr(pauli(n,j)'*M))
        end
    end
end

#Check that coefficient distribution agrees with prediction

#Poisson statistics
x = 0:n+1
poisson(x) = 1/exp(1)*1/factorial(x)
dist_X
hist = sort(countmap(round.(dist_X, digits = 5)))
bar(collect(keys(hist)), collect(values(hist)./length(dist_X)))
scatter!(x./sqrt(D), poisson.(x), markersize = 5)
ylims!(0, .4)

#Bessel statistics
x = -n:n
modified_bessel(x) = 1/exp(1) * sum([1/(gamma(n+1)gamma(n+abs(x)+1))*(1/2)^(2n+abs(x)) for n in 0:10])
hist = sort(countmap(round.(dist_XZ, digits = 5)))
bar(collect(keys(hist)), collect(values(hist)./length(dist_XZ)))
scatter!(x./sqrt(D), modified_bessel.(x), markersize = 5)

#Try going the other way: sample based on these statistics and bin singular values
function sample_pauli_perm(D)
    n = round(Int, log2(D)) #total number of qubits
    w_poisson = Weights(poisson.(0:n))
    w_bessel = Weights(modified_bessel.(-n:n))
    mat = zeros(D,D)
    for i in 1:D
        for j in 1:D
            if i ≤ sqrt(D) && j ≤ sqrt(D)
                mat[i,j] = sample((0:n)./sqrt(D), w_poisson)
            else
                mat[i,j] = sample((-n:n)./sqrt(D), w_bessel)
            end
        end
    end
    mat
end

binomial(n, p, N) = gamma(N+1)/(gamma(N-n+1)*gamma(n+1))*p^n*(1-p)^(N-n)

#try adding in dependence of the variance
function sample_pauli_perm_correlated(D)
    n = round(Int, log2(D)) #total number of qubits
    x_strings = round(Int, sqrt(D))
    w_poisson = Weights(poisson.(0:n))
    mat = zeros(x_strings, x_strings, x_strings, x_strings)
    for i in 1:x_strings
        for j in 1:x_strings
            r = sample((0:n), w_poisson)
            mat[i,1,j, 1] = r
            if r ≠ 0
                w_binom = Weights([binomial(k, 1/2, r) for k in 0:r])
            end
            for k in 1:x_strings
                for l in 1:x_strings
                    if r ≠ 0
                        mat[i,k,j,l] = sample((0:r), w_binom)/sqrt(D)
                    end
                end
            end

       end
    end
    reshape(mat, (D,D))
end


#Resampling does not produce the same statistics as a random permutation!
dist = []
for i in ProgressBar(1:100)
    M = sample_pauli_perm(2^8)
    U,Λ,V = svd(M)
    map(x->push!(dist, x), Λ)
end

histogram(dist) #WD statistics are the limit of free convolution, like the free analogue of a Gaussian

             
#still semi-circle distributed!
dist = []
for i in ProgressBar(1:50)
    M = sample_pauli_perm(2^10)
    U,Λ,V = svd(M)
    map(x->push!(dist, x), Λ)
end

histogram(dist, normalize = true)
xlims!(0, 2)

#Trying to model atoms with Gaussian / Poisson...
dist = []
dist2 = []
for i in ProgressBar(1:10)
    M = rand_perm_mat(2^12)
    d = ES(M, 2^6)
    d = round.(d, digits = 5)
    map(x->push!(dist, x), d)
    map(x->push!(dist2, x), diag(M'*M))
end

poisson(x, q) = exp(-q)*q^x/gamma(x+1)
xax = LinRange(0,7, 100)

emp_bins = -.5:1:7.5
cnmp = countmap(dist)
cnmp2 = countmap(dist2)
diag_dist = [cnmp2[i]/sum(values(cnmp2)) for i in 0:7]

#attempt #1: histogram binning at integer
histogram(dist, bins = emp_bins, normalize = true, label = "ES histogram", linetype = :stephist, linewidth = 3, dpi = 250)
bar!(-.5:1:6.5, diag_dist, linewidth = 3, linetype = :step, label = L"diags $M^\dagger M$")
scatter!(0:7, x->poisson(x, 1), label = "Pois(1)")
title!(L"Histogram of ES vs diagonal values of $M^\dagger M$")
xlabel!(L"\lambda^2")
savefig("figs/ESvsDiags.png")

#Attempt 2: isolate the point spectrum
point_spec = [float(cnmp[i]) for i in 0:5]
point_spec ./= sum(point_spec)
bar(-.5:4.5, point_spec, linetype = :step, label = "ES integral point spectrum", linewidth = 3, dpi = 250)
bar!(-.5:1:6.5, diag_dist, linewidth = 3, linetype = :step, label = L"diags $M^\dagger M$")
scatter!(0:6, x->poisson(x, 1), label = "Pois(1)")
title!(L"ES integral point spectrum vs diagonal values of $M^\dagger M$")
xlabel!(L"\lambda^2")
savefig("figs/ESPSvsDiags.png")

#Attemp 3: neglect zeros
point_spec = [float(cnmp[i]) for i in 1:5]
point_spec ./= sum(point_spec)
diag_dist = [cnmp2[i]/sum(values(cnmp2)) for i in 1:7]
diag_dist ./= sum(diag_dist)
pois_dist = [poisson(x, 1) for x in 1:5]
pois_dist ./= sum(pois_dist)
pois_dist2 = [poisson(x, 1-exp(-1)) for x in 1:5]
pois_dist2 ./= sum(pois_dist2)
bar(.5:4.5, point_spec, linetype = :step, label = "ES point spectrum (renormalized)", linewidth = 3, dpi = 250)
bar!(.5:1:6.5, diag_dist, linewidth = 3, linetype = :step, label = L"diags $M^\dagger M$")
scatter!(1:5, pois_dist, label = "Pois(1) (renorm.)")
scatter!(1:5, pois_dist2, label = "Pois(1-1/e) (renorm)")
title!("Non-zero ES")
xlabel!(L"\lambda^2")
savefig("figs/ESVsPSNonzero.png")


#Now try making a cut at different places. How does the spectrum change?
dists = []
N = 10
Cs = 1:5

for c in Cs
    Da = 2^c
    Db = 2^(N-c)
    dist = []
    for s in ProgressBar(1:2^(N-c))
        rperm = rand_perm_mat(2^N)
        U, D, V = tsvd(TensorMap(rperm, ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db), (1,3), (2,4))
        map(x->push!(dist, x), diag(D.data))
    end
    push!(dists, dist)
end

p = plot(dpi = 250)
title!("Permutation entanglement spectrum by cut length")
for i in 1:5
    histogram!(dists[i], normalize = true, label = "$(i)")
    vline!([2^(10/2-i)], label = "")
end
xlabel!("λ")
ylabel!("PDF")
display(p)
savefig("figs/Cut_size_vs_ES.png")

p = plot(dpi = 250)
title!("Permutation entanglement spectrum by cut length")
for i in 1:5
    cm = countmap(dists[i])
    bar!(collect(keys(cm)), collect(values(cm))/sum(values(cm)), label = "", linewidth = .5)
    vline!([2^(10/2-i)], label = "")
end
xlabel!("λ")
ylabel!("frac")
display(p)
savefig("figs/Cut_size_vs_ES_bar.png")


cm = countmap(dists[5])
cm2 = countmap(round.(dist, digits = 10))
bar(collect(keys(cm2)), collect(values(cm2))./sum(values(cm2)))
bar!(collect(keys(cm)), collect(values(cm))./sum(values(cm))*15)

Da = 2
Db = 2


fixedpts(p) = sum(length.(cycle_notation(p)) .== 1)

dist = []
for i in ProgressBar(1:8000)
    P1 = rand_perm_mat(2^8)
    P2 = rand_perm_mat(2^8)
    M = [tr(P1^2)+((-1)^i+(-1)^j)*tr(P1*P2)+(-1)^(i+j)*tr(P2^2) for i in 1:2, j in 1:2]
    e = svdvals(M)
    push!(dist, e[1])
    push!(dist, e[2])
end

histogram(dist)

Da = 2
Db = 2
dist = []
for s in ProgressBar(1:8000)
    rperm = rand_perm_mat(Da*Db)
    U, D, V = tsvd(TensorMap(rperm, ℂ^Da*ℂ^Db, ℂ^Da*ℂ^Db), (1,3), (2,4))
    map(x->push!(dist, x), diag(D.data))
end

bar(dist, noramlize = true)


histogram(dist)



cycle_notation(P)
fixedpts(P)

P1 = rand_perm_mat(2)
P2 = rand_perm_mat(2)
M = [tr(P1^2)+((-1)^i+(-1)^j)*tr(P1*P2)+(-1)^(i+j)*tr(P2^2) for i in 1:2, j in 1:2]
M



#Next, trying to see if obsvious patterns emerge in X^†X

#For Bernoulli...

D = 2^6
X = convert.(Float64,rand(D, D) .< 1/D)
heatmap(X'*X)

# Can we sample X^\dag X randomly?
N = 2^12

w1 = Weights([exp(-1)/factorial(x) for x in 0:5])
w2 = Weights([exp(-1/N)*1/N^x/factorial(x) for x in 0:2])

M = zeros(Float64, N,N)
for i in 1:N
    for j in 1:i
        if i != j
            M[i,j] = sample(0:2, w2)
            M[j,i] = M[i,j]
        else
            M[i,i] = sample(0:5, w1)
        end
    end
end

M
p = heatmap(M)
display(p)

vals, vecs = eigen(M)
mp = countmap(vals)
histogram(vals) #not positive, and produces a very different distribution!
