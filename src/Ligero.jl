module Ligero

using BinaryFields, MerkleTree, BinaryReedSolomon

is_pow_two(x) = x & (x - 1) == 0
next_pow_two(x) = 1 << (1 + floor(Int, log2(x - 1)))

# function evaluate_lagrange_basis(rs::Vector{T}) where T <: BinaryElem
#     one_elem = T(one(T))
#     current_layer = [one_elem + rs[1], rs[1]]
#     len = 2
#     for i in 2:length(rs)
#         next_layer_size = 2 * len
#         next_layer = Vector{T}(undef, next_layer_size)

#         ri_p_one = one_elem + rs[i]
#         for j in 1:len
#             next_layer[2*j - 1] = current_layer[j] * ri_p_one
#             next_layer[2*j]   = current_layer[j] * rs[i]
#         end

#         current_layer = next_layer
#         len *= 2
#     end

#     return current_layer
# end

# struct Message{T <: BinaryElem}
#     data::Matrix{T}
# end

struct LigeroProofProperties
    log2_prob::Int
    inv_rate::Int
    base_field_size::Int
    ext_field_size::Int
end

n_queries(; log2_prob=100, inv_rate=4) = ceil(Int, -log2_prob/log2((1+1/inv_rate)/2))

# Assumes `N` is a power of 2
function opt_dims(N, prop)
    S = n_queries(; prop.log2_prob, prop.inv_rate)
    n = next_pow_two(sqrt(prop.ext_field_size*N/(prop.base_field_size*S)))

    return div(N, n), n
end

struct LigeroProverCommitment{T}
    mat::Matrix{T}
    tree::Vector{Vector{String}}
end

struct LigeroVerifierCommitment{T}
    root::String
end

function commit(poly::Vector{T}; properties=nothing) where T <: BinaryElem
    @assert length(poly) > 0 "Polynomial must have at least one element"
    @assert is_pow_two(length(poly)) "Polynomial length must be a power of 2"

    if isnothing(properties)
        log2length_poly = round(Int, log2(length(poly)))
        if log2length_poly <= 14
            properties = LigeroProofProperties(100, 4, 16, 128)
        else
            properties = LigeroProofProperties(100, 4, 32, 128)
        end
    end

    m, n = opt_dims(length(poly), properties)
    poly_mat = reshape(poly, m, n)
    # Do the encoding here

    leaves = [BinaryFields.binary_val.(row) for row in eachrow(poly_mat)]
    tree = build_merkle_tree(leaves)

    return LigeroProverCommitment(poly_mat, tree), LigeroVerifierCommitment(tree[end][1])
end

# struct InterleavedCode{T <: BinaryElem}
#     code::Matrix{T}
# end
# export InterleavedCode

# function InterleavedCode(data::Matrix{T}, column_length::Int) where T <: BinaryElem
#     if size(data, 1) != column_length
#         error("Invalid code: Expected each column to have length $column_length, but got $(size(data, 1)).")
#     end
#     InterleavedCode{T}(data)
# end

# function encode(msg::Message{T}, rs::ReedSolomonEncoding{T}) where {T <: BinaryElem}
#     @assert size(msg.data, 1) == message_length(rs) "Each message column must have length $(message_length(rs)), got $(size(msg.data, 1))"
#     encoded_columns = [encode(rs, Vector(col)) for col in eachcol(msg.data)]
#     encoded_matrix = hcat(encoded_columns...)
#     return InterleavedCode(encoded_matrix)
# end

# function commit(code::InterleavedCode{T}) where {T <: BinaryElem}
#     leaves = [map(BinaryFields.binary_val, Vector(row)) for row in eachrow(code.code)]
#     return build_merkle_tree(leaves)
# end

# # TODO: convert from GF_2^16 to GF_2^128
# function prove(msg::Message{T}, gr::Vector{T}) where {T <: BinaryElem}
#     @assert size(msg.data, 2) == length(gr) "Number of columns must match length of gr"
#     return msg.data * gr
# end

# export encode, commit, prove, evaluate_lagrange_basis

end
