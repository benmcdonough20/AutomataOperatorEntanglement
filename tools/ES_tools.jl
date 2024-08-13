using StatsBase
using Random

function ptranspose!(M, hsize)
	for i in 0:hsize^2-1
	    for j in 0:hsize^2-1
		in2 = i % hsize
		in1 = i ÷ hsize
		out2 = j % hsize
		out1 = j ÷ hsize
		if in2 < out1
		    tmp = M[i+1,j+1]
		    M[i+1,j+1] = M[in1*hsize+out1+1, in2*hsize + out2+1]
		    M[in1*hsize+out1+1, in2*hsize+out2+1] = tmp
		end
	    end
	end
end

function ptranspose!(M)
	hsize = Int(sqrt(size(M, 1)))
	for i in 0:hsize^2-1
	    for j in 0:hsize^2-1
		in2 = i % hsize
		in1 = i ÷ hsize
		out2 = j % hsize
		out1 = j ÷ hsize
		if in2 < out1
		    tmp = M[i+1,j+1]
		    M[i+1,j+1] = M[in1*hsize+out1+1, in2*hsize + out2+1]
		    M[in1*hsize+out1+1, in2*hsize+out2+1] = tmp
		end
	    end
	end
end

ES = (op, D) -> (
    ptranspose!(op, D);
    svdvals(op).^2
)

ES2 = (op) -> (
    ptranspose!(op);
    svdvals(op)
)

S(E) = sum(map(x->-x*log(x), E[E.>0]))

function merge_countmaps!(cm1, cm2)
    for (key, val) in cm2
        if !haskey(cm1, key)
            cm1[key] = 0
        end
        cm1[key] += val
    end
end

function add_countmap!(dict, key, val)
    if !haskey(dict, key)
        dict[key] = 0
    end
    dict[key] += val
end

function rand_PHP(N, nH)
	H = 1/sqrt(2) * [1 1; 1 -1]
	rand_locs = sample(1:N, nH, replace = false)
	Hf = insert_op_at_locs(H, N, rand_locs)
	Base.permutecols!!(Base.permutecols!!(Hf, randperm(2^N))', randperm(2^N))'
end

function spacings(dist; adj = false, prec = 0)
	if prec != 0	
		dist = round.(dist, digits = prec)
	end
	d = reverse(sort(dist));
	arr = [(d[i+2]-d[i+1])/(d[i+1]-d[i]) for i in 1:length(dist)-2]
	if adj
		arr = [min(r, 1/r) for r in arr]
	end
	arr
end