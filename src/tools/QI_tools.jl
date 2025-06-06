function insert_op_at_locs(op, num_qubits, locs)
    return kron([(i in locs) ? op : Matrix(I, 2, 2) for i in 1:num_qubits]...)
end

function level_spacings(dist)
    d=sort(dist)
    ret = [(d[i]-d[i-1])/(d[i-1]-d[i-2]) for i in 3:length(d)]
    ret[abs.(ret) .<= 25]
end