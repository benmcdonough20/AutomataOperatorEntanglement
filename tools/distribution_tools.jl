semicircle(x) = 1/(2π*sqrt(x))*sqrt(4-x)
semicircle_CDF(x) = x ≥ 4 ? 1 : (x ≤ 0 ? 0 : 1/2π*(2π+x*sqrt(4/x-1)-4*atan(sqrt(4/x-1))))

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
    step = (stop-start)/num
    [sum(start+(i-1)*step .<= dist .< start + i*step) for i in 1:num-1] ./ (length(dist) * step)
end