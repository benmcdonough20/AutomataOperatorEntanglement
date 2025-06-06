#represent a Pauli operator on n qubits as X string followed by Z string
struct pauli_string
    n::Int
    X::Vector{Bool}
    Z::Vector{Bool}
end

bin_to_int(bin) = sum([b*2^(i-1) for (i,b) in enumerate(reverse(bin))])
int_to_bin(num, dim) = [(num >> i)%2 for i in 0:dim-1]
pauli_str_to_int(pauli_str) = bin_to_int(vcat(pauli_str.Z, pauli_str.X))
int_to_pauli_str(num, dim) = (x=int_to_bin(num, dim*2);pauli_string(dim, x[1:dim],x[dim+1:end]))

const X_op = [0 1.0; 1.0 0]
const Z_op = [1.0 0 ; 0 -1.0]

#building block Pauli operators 00 01 10 11
function primitive_pauli(x::Bool, z::Bool, real = false)
    1/sqrt(2)*(real ? 1 : 1im)^(x*z)*X_op^x*Z_op^z
end

##Create a Pauli operator from a string
function pauli(pauli_str::pauli_string, real = false)
    if pauli_str.n == 1
        primitive_pauli(first(pauli_str.X), first(pauli_str.Z), real)
    else
        kron([primitive_pauli(x,y, real) for (x,y) in zip(pauli_str.X, pauli_str.Z)]...)
    end
end

#create a Pauli oeprator from an index
function pauli(nqubits::Int, idx::Int, real = false)
    @assert idx < 4^nqubits
    return pauli(int_to_pauli_str(idx, nqubits), real)
end

#vectorize a matrix in the Pauli basis
function pauli_basis(M, nqubits, real = false)
    @assert first(size(M)) == last(size(M)) == 2^nqubits
    return [tr(pauli(nqubits, i, real)'*M) for i in 0:4^nqubits-1]
end