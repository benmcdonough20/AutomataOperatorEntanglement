using Graphs
using GraphPlot
using ProgressBars
using SpecialFunctions
using Polylogarithms
include("tools/imports.jl")


include("tools/plotting_defaults.jl")


#Note: lots of variability with more nodes and increased samples...
#resolution: probability that node belongs to component of size k is not the same as the 
#probability that any given component will have size k
N = 2^12
samples = 100
dist = countmap([])

#Let's check that this analysis also works for permutations
X = kron([0 1; 1 0], Matrix(I, Int(N/2), Int(N/2)))
for i in ProgressBar(1:samples)
    M = rand_perm_mat(N)
    #M = M*X*M'
    ptranspose!(M)
    g = SimpleDiGraph(M)
    for v in 1:N
        for cmp in weakly_connected_components(g)
            if v in cmp
                add_countmap!(dist, length(cmp), 1)
            end
        end
    end
end

x = LinRange(1,10, 100)
bar(collect(keys(dist)), collect(values(dist))./sum(values(dist)), label = "emp", dpi = 250, linewidth = 0, color = :blue)
title!("Vertex component sizes in random automata")
psmall = sum([exp(-2k)*(2k)^(k-1)/gamma(k+1) for k in 1:10])
1-psmall
vline!([real(1-psmall)*N], label = L"1-p_{small}")
xlabel!("|C(v)|")
ylabel!("P(|C(v)|)")
savefig("figs/paper/vertex_cmp_sizes_permutation.svg")

xlims!(1, 10)
ylims!(0,.15)
c = 2
plot!(x, x->exp(-c*x)*(c*x)^(x-1)/gamma(x+1), label = "Poisson branching process", linewidth = 3, )
savefig("figs/vertex_cmp_sizes_permutation_small.svg")

#specific example to illustrate appearance of one large connected component
N = 2^10
M = rand_perm_mat(N)
dim = round(Int, sqrt(N))
T = TensorMap(M, ℂ^dim * ℂ^dim, ℂ^dim*ℂ^dim)
@tensor Tpt[a,b,c,d] := T[a,c,d,b]
M .= reshape(Tpt.data, (N,N))
g = SimpleDiGraph(M)
weakly_connected_components(g)
ctmap = countmap(length.(weakly_connected_components(g)))
bar(collect)



#Hypothesis: The number of components with only one node should be the umber of nodes with no incoming or outgoing edges
#This will happen with probability 1/e^2
#This agrees with the Poisson branching process estimation for the case c = 2

# self-connected elements: N*1/eN - expected number of 1's - 1/e, so expected fraction on 1's is 1/eN

#Poisson branching process
#Page 181 of probabilistic methods


#Joint distribution of large and small pieces! Check Erdos-Renyi phase transition, methods ch 11

#What if we look at the distribution of cycles of length n in X^†X?

#Could we find the probability distribution over component sizes for small k?

#Wolfram gives sum of the above as PolyLog(3/2, 2/e),
#so the probability of belonging to a "small" component is around 1/2

N = 2^10
samples = 100
small_sizes = 0

for i in ProgressBar(1:samples)
    #M = rand(N,N) .< 1/N
    M = rand_perm_mat(N)
    ptranspose!(M)
    g = SimpleDiGraph(M)
    for v in 1:N
        for cmp in weakly_connected_components(g)
            if v in cmp
                if length(cmp) < 20
                    small_sizes += 1
                end
            end
        end
    end
end

small_sizes / (samples*N)
sum([exp(-c*x)*(c*x)^(x-1)/gamma(x+1) for x in 1:100])
sum([1/sqrt(2π)*exp(-x)*2^(x-1)*x^(-3/2) for x in 1:100])
# 1/2√2π PolyLog(3/2, 1/e)

dist = countmap([])
N = 2^13
samples = 10

for i in ProgressBar(1:samples)
    M = rand(N,N) .< 1/N
    g = SimpleDiGraph(M)
    add_countmap!(dist, sum(length.(weakly_connected_components(g)) .> 100), 1)
end

bar(collect(keys(dist)), collect(values(dist)))
dist
#Almost always only one large component, I wonder if this holds in general...
#Yay this considerably simplifies the analysis. This makes it easy to predict the distribution of components on the smaller end, and also on the bigger end (although need an analytic reason for only 1!)
#polylog function... also arises in fermi-dirac distribution

#So it follows that distribution of component size is <1/x>

#log gammas
#exponentiate log gamma

N = 2^9
dist =[]
samples = 1000
#Now look specifically at spectrum of connected component
for i in ProgressBar(1:samples)
    M = convert.(Float64, rand(N,N) .< 1/N)
    g = SimpleDiGraph(M)
    cmps = weakly_connected_components(g)
    big_vs = []
    for n in 1:N
        for cmp in cmps 
            if length(cmp) > .5*N && n in cmp
                push!(big_vs, n)
            end
        end
    end
    M = M[big_vs,big_vs]
    U,D,V = svd(M)
    #map(x->push!(dist, x), Diagonal(D)[round.(Diagonal(D),digits = 10) .!= 0 ])
    map(x->push!(dist, x), D)
end

mp = countmap(round.(dist, digits = 5))
empbins = LinRange(0,3,83)
histogram(dist, normalize = true, bins = empbins, label = "", dpi = 300)
xlabel!(L"$\lambda$")
ylabel!("PDF")
savefig("figs/big_comp_hist.png")

bar(collect(keys(mp)), collect(values(mp))./sum(values(mp)), label = "")

xlabel!(L"$\lambda$")
ylabel!("PDF")
savefig("figs/big_comp_hist.png")
vline!([sqrt(1), sqrt(2), sqrt(3), 1.5])
sort(collect(keys(mp)), lt = (x,y)->mp[x] > mp[y])

N = 2^7
M = convert.(Float64, rand(N,N) .< 1/N)
g = SimpleDiGraph(M)
gplot(g)
cmps = weakly_connected_components(g)
big_vs = []
for n in 1:N
    for cmp in cmps 
        if length(cmp) > .5*N && n in cmp
            push!(big_vs, n)
        end
    end
end

M = M[big_vs,big_vs]
gplot(SimpleDiGraph(M))

heatmap(M'*M)

sum([(M'*M)[20,i] for i in 1:size(M)[1] if i ≠ 20])

V,E = eigen(M'*M);
V = round.(V, digits = 10)
Ones = hcat([E[:,i] for i in 1:size(M)[1] if V[i] ≈  1]...)
#2.6180339887 shows up in big graph and small graph... why?
V
heatmap(Ones)
Ones[:,1]' * M'*M * Ones[:,1]

using Colors
c2 = colorant"black"
c1 = colorant"red"
col = x->[weighted_color_mean((y-minimum(x))/(maximum(x)-minimum(x)), c1, c2) for y in x]
heatmap(Ones)
gplot(SimpleDiGraph(M'), nodefillc = col(round.(abs.(Ones[:,1]), digits = 5) .> 0))#, nodelabel = 1:size(M)[1], nodelabeldist = 1.5)
gplot(SimpleDiGraph(M), nodefillc = col(Ones[:,1]))
M
Ones
Ones[:,1]
norm(M'*M * Ones[:,1]-Ones[:, 1])
V

M
barplt(V)
maximum(svdvals(M))
histogram(svdvals(M))
barplt(svdvals(M))

M*E[2,:] == E[2,:]
V[V .!= 0]
E[1,:][E[1,:] .!= 0]

using Colors


D
M

gp = SimpleDiGraph(M')
gplot(gp)

eigen(M)

1/sqrt(2)
1/sqrt(3)
sqrt(2/5)

barplt(dist) = (
    cm = countmap(dist);
    bar(collect(keys(cm)), collect(values(cm))/sum(values(cm)))
)

N = 500
M = convert.(Float64, rand(N,N) .< 1/N)
g = SimpleDiGraph(M)
cmps = weakly_connected_components(g)
big_vs = []
for n in 1:N
    for cmp in cmps 
        if length(cmp) > .5*N && n in cmp
            push!(big_vs, n)
        end
    end
end

M = M[big_vs,big_vs]
D = svdvals(M)
#map(x->push!(dist, x), Diagonal(D)[round.(Diagonal(D),digits = 10) .!= 0 ])
dist = []
map(x->push!(dist, x), D)

barplt(dist)