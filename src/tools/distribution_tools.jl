semicircle(x) = 1/(2π*sqrt(x))*sqrt(4-x)
semicircle_CDF(x) = x ≥ 4 ? 1 : (x ≤ 0 ? 0 : 1/2π*(2π+x*sqrt(4/x-1)-4*atan(sqrt(4/x-1))))
wigner_PDF(x) = 0 <= x <= 2 ? 1/π * sqrt(4-x^2) : 0 #Integration convention: use right sums!
wigner_CDF(x) = 0 <= x <= 2 ? 1/π * (x/2*sqrt(4-x^2) + 2 * asin(x/2)) : (x < 0 ? 0 : 1)
WD_PDF(x) = 27/8 * (x+x^2)/(1+x+x^2)^(5/2)
WD_PDF_adj(x) = WD_PDF(x) + B/(1+x)^2*(1/(x+1/x)-A*1/(x+1/x)^2)
WD_CDF(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2))

B = 0.233378
A = 2*(π-2)/(4-π)
δWD(x) = 1/4*B*(-2*(1+A)+(2+A)/(1+x)+A/(1+x^2)+(2+A)*atan(x))
WD_CDF_adj(x) = 1/4*(2+(1+2x)*(x^2+x-2)/(1+x+x^2)^(3/2)) + δWD(x)

struct distribution
    range
    PDF
    CDF
end

struct bins
    range
    w
end

semicircle_distribution = distribution((0, 4), semicircle, semicircle_CDF)

function binlist(bins)
    collect(first(bins.range):bins.w:last(bins.range))
end

function axes(bins)
    binlist(bins)[1:end-1] .+ bins.w/2
end

function discretize(dist::distribution, bin::bins)
    bin_list = binlist(bin)
    [(dist.CDF(bin_list[i+1]) - dist.CDF(bin_list[i]))/bin.w for i in 1:(length(bin_list)-1)]
end

function hist(dist::Vector, bin)
    bin_list = binlist(bin)
    if minimum(dist) < first(bin_list) || maximum(dist) > last(bin_list)
        throw(DomainError((minimum(dist), maximum(dist)), "Increase bin range!"))
    end
    [sum((x->bin_list[i]<x≤bin_list[i+1]).(dist))/(length(dist)*bin.w) for i in 1:length(bin_list)-1]
end

function KLD(Q, actual_dist::distribution, bin)
    P = discretize(actual_dist, bin)
    if any(Q[P .!= 0] .== 0)
        throw(DomainError(0, "measures cannot be mutually singular"))
    end
    sum([p==0 ? 0 : p*(log(p)-log(q)) for (p,q) in zip(P, Q)])*bin.w
end

function Hell(Q, actual_dist::distribution, bin)
    P = discretize(actual_dist, bin)
    sum([(p-q)^2 for (p,q) in zip(P, Q)])*bin.w
end

function Wass(Q, actual_dist::distribution, bin)
    cdf2 = [sum(Q[1:i]) for i in 1:length(Q)]
    sum(abs.([cdf2 .- actual_dist.CDF(i * bin.w + range(0)) for i in 0:length(Q)-1]))*bin.w
end

function bin(dist, start, stop, num)
    x = LinRange(start, stop, num)
    [sum(x[i] .<= dist .< x[i+1]) for i in 1:length(x)-1]/length(dist) * (num-1) / stop
end