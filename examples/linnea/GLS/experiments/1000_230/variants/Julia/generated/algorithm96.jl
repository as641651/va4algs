using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm96(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm96(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.76e+08 FLOPs
    # X: ml0, full, M: ml1, full, y: ml2, full
    ml3 = Array{Float64}(undef, 230, 1000)
    # tmp56 = X^T
    transpose!(ml3, ml0)

    # X: ml0, full, M: ml1, full, y: ml2, full, tmp56: ml3, full
    # (L2 L2^T) = M
    LAPACK.potrf!('L', ml1)

    # X: ml0, full, y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular
    # tmp12 = (L2^-1 X)
    trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0)

    # y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular, tmp12: ml0, full
    ml4 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    syrk!('L', 'T', 1.0, ml0, 0.0, ml4)

    # y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular, tmp14: ml4, symmetric_lower_triangular
    # tmp68 = (L2^-1 y)
    trsv!('L', 'N', 'N', ml1, ml2)

    # tmp56: ml3, full, L2: ml1, lower_triangular, tmp14: ml4, symmetric_lower_triangular, tmp68: ml2, full
    # (L15 L15^T) = tmp14
    LAPACK.potrf!('L', ml4)

    # tmp56: ml3, full, L2: ml1, lower_triangular, tmp68: ml2, full, L15: ml4, lower_triangular
    # tmp71 = (L2^-T tmp68)
    trsv!('L', 'T', 'N', ml1, ml2)

    # tmp56: ml3, full, L15: ml4, lower_triangular, tmp71: ml2, full
    # tmp219 = (L15^-1 tmp56)
    trsm!('L', 'L', 'N', 'N', 1.0, ml4, ml3)

    # L15: ml4, lower_triangular, tmp71: ml2, full, tmp219: ml3, full
    ml5 = Array{Float64}(undef, 230)
    # tmp23 = (tmp219 tmp71)
    gemv!('N', 1.0, ml3, ml2, 0.0, ml5)

    # L15: ml4, lower_triangular, tmp23: ml5, full
    # tmp24 = (L15^-T tmp23)
    trsv!('L', 'T', 'N', ml4, ml5)

    # tmp24: ml5, full
    # b = tmp24
    return (ml5)
end