using StatsBase
using Ligero
using BinaryFields, BinaryReedSolomon

poly = rand(BinaryElem16, 2^24)

@info "Constructing commitment"
comm = commit(poly; verbose=true, parallel=false)
verifier_comm = verifier_commitment(comm)
@info "Verifier commitment size: $(Base.format_bytes(sizeof(verifier_comm)))"

# For now, assume interactive randomness
n_rows, n_cols = size(matrix(comm))
gr = rand(BinaryElem128, n_cols)

# 148 is the number of necessary queries for soundness 2^(-100)
# for code with rate 1/4
S_sorted = sort(sample(1:n_rows, 148, replace=false))

@info "Constructing proof"
@time proof = prove(comm, gr, S_sorted)
proof_size = sizeof(proof)
@info "Proof size: $(Base.format_bytes(proof_size))"

@info "Verifying proof"
verified = verify(proof, verifier_comm, S_sorted, gr)
@show verified;