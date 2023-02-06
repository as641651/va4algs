using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm83(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm83(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.75e+08 FLOPs
    # X: ml0, full, M: ml1, full, y: ml2, full
    # (L2 L2^T) = M
    LAPACK.potrf!('L', ml1)

    # X: ml0, full, y: ml2, full, L2: ml1, lower_triangular
    # tmp68 = (L2^-1 y)
    trsv!('L', 'N', 'N', ml1, ml2)

    # X: ml0, full, L2: ml1, lower_triangular, tmp68: ml2, full
    # tmp12 = (L2^-1 X)
    trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0)

    # tmp68: ml2, full, tmp12: ml0, full
    ml3 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    gemm!('T', 'N', 1.0, ml0, ml0, 0.0, ml3)

    # tmp68: ml2, full, tmp12: ml0, full, tmp14: ml3, full
    ml4 = Array{Float64}(undef, 230)
    # tmp21 = (tmp12^T tmp68)
    gemv!('T', 1.0, ml0, ml2, 0.0, ml4)

    # tmp14: ml3, full, tmp21: ml4, full
    # (L15 L15^T) = tmp14
    LAPACK.potrf!('L', ml3)

    # tmp21: ml4, full, L15: ml3, lower_triangular
    # tmp23 = (L15^-1 tmp21)
    trsv!('L', 'N', 'N', ml3, ml4)

    # L15: ml3, lower_triangular, tmp23: ml4, full
    # tmp24 = (L15^-T tmp23)
    trsv!('L', 'T', 'N', ml3, ml4)

    # tmp24: ml4, full
    # b = tmp24
    return (ml4)
end