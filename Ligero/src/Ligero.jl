module Ligero

import Base: sizeof

using Base.Threads, ThreadTools
using BinaryFields, BatchedMerkleTree, BinaryReedSolomon

export LigeroProofProperties, LigeroCommitment, LigeroVerifierCommitment, LigeroProofProperties
export commit, matrix, verifier_commitment, prove, verify

is_pow_two(x) = x & (x - 1) == 0
next_pow_two(x) = 1 << (1 + floor(Int, log2(x - 1)))

struct LigeroProofProperties
    log2_prob::Int
    inv_rate::Int
    base_field_size::Int
    ext_field_size::Int
end

n_queries(; log2_prob=100, inv_rate=4) = ceil(Int, -log2_prob/log2((1+1/inv_rate)/2))

# Assumes `N` is a power of 2
function opt_dims(N, prop::LigeroProofProperties)
    S = n_queries(; prop.log2_prob, prop.inv_rate)
    n = next_pow_two(sqrt(prop.ext_field_size*N/(prop.base_field_size*S)))

    return div(N, n), n
end

struct LigeroCommitment{T}
    mat::Matrix{T}
    tree::CompleteMerkleTree
    rs::BinaryReedSolomon.ReedSolomonEncoding{T}
end

matrix(c::LigeroCommitment) = c.mat

struct LigeroVerifierCommitment
    root::MerkleRoot
end

sizeof(x::LigeroVerifierCommitment) = sizeof(x.root)

verifier_commitment(com) = LigeroVerifierCommitment(get_root(com.tree))

# XXX: Make this allocate way less (a lot of GC pressure)
function encode_cols(poly_mat, rs; parallel=true)
    if parallel
        encoded_columns = tmap(c -> encode(rs, c), eachcol(poly_mat))
    else
        encoded_columns = map(c -> encode(rs, c), eachcol(poly_mat))
    end
    return hcat(encoded_columns...)
end

"""
    commit(poly; properties=nothing)

Commit to a polynomial `poly`, represented as an AbstracVector of `BinaryElem16` or
`BinaryElem32`.
"""
function commit(poly; properties=nothing, verbose=false, parallel=true)
    @assert length(poly) > 0 "Polynomial must have at least one element"
    @assert is_pow_two(length(poly)) "Polynomial length must be a power of 2"

    T_poly = eltype(poly)

    if isnothing(properties)
        properties = LigeroProofProperties(100, 4, bitsize(T_poly), 128)
    end

    m, n = opt_dims(length(poly), properties)
    min_field_size = round(Int, log2(m))+round(Int, log2(properties.inv_rate))
    if min_field_size > properties.base_field_size
        error("Base field size $(properties.base_field_size) is too small for polynomial of length $(length(poly)) with inv_rate $(properties.inv_rate). Minimum size is $min_field_size.")
    end
    rs = reed_solomon(T_poly, m, m*properties.inv_rate)

    poly_mat = reshape(poly, m, n)
    if verbose
        @info "Encoding"
        @time mat = encode_cols(poly_mat, rs; parallel)
    else
        mat = encode_cols(poly_mat, rs)
    end

    leaves = eachrow(mat)
    
    if verbose
        @info "Building Merkle tree"
        @time tree = build_merkle_tree(leaves)
    else
        tree = build_merkle_tree(leaves)
    end

    return LigeroCommitment(mat, tree, rs)
end

struct LigeroProof{T, U} 
    y_r::Vector{U}
    merkle_openings::BatchedMerkleProof
    rows::Matrix{T}
    rs::BinaryReedSolomon.ReedSolomonEncoding{T}
end

sizeof(x::LigeroProof) = sizeof(x.y_r) + sizeof(x.merkle_openings) + sizeof(x.rows)

function prove(comm::LigeroCommitment{T}, gr, S_sorted) where {T <: BinaryElem}
    n_rows = message_length(comm.rs)

    openings = MerkleTree.prove(comm.tree, S_sorted)
    y_r = @views comm.mat[1:n_rows, :] * gr

    return LigeroProof(y_r, openings, comm.mat[S_sorted, :], comm.rs)
end

function verify(proof::LigeroProof, com::LigeroVerifierCommitment, S_sorted, gr)
    t_depth = log_block_length(proof.rs)

    # Merkle inclusion proof verification
    # XXX: maybe improve output, specify which row proof failed
    if !MerkleTree.verify(com.root, proof.merkle_openings;
            depth = t_depth,
            leaves = eachrow(proof.rows),
            leaf_indices = S_sorted
        )
        @warn "Merkle inclusion proof failed to verify"
        return false
    end

    # Looks very innocent! Consider the fields being used here
    enc_y = encode(proof.rs, proof.y_r)
    expected_output = @views proof.rows * gr

    # Check if the output of the encoding matches the expected output
    if enc_y[S_sorted] != expected_output
        @warn "Encoding verification failed"
        return false
    end

    return true
end


end
