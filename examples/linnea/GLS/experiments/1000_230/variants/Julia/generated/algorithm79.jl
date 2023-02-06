using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm79(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm79(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.75e+08 FLOPs
    # X: ml0, full, M: ml1, full, y: ml2, full
    # (L2 L2^T) = M
    LAPACK.potrf!('L', ml1)

    # X: ml0, full, y: ml2, full, L2: ml1, lower_triangular
    # tmp12 = (L2^-1 X)
    trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0)

    # y: ml2, full, L2: ml1, lower_triangular, tmp12: ml0, full
    # tmp68 = (L2^-1 y)
    trsv!('L', 'N', 'N', ml1, ml2)

    # tmp12: ml0, full, tmp68: ml2, full
    ml3 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    syrk!('L', 'T', 1.0, ml0, 0.0, ml3)

    # tmp12: ml0, full, tmp68: ml2, full, tmp14: ml3, symmetric_lower_triangular
    # (L15 L15^T) = tmp14
    LAPACK.potrf!('L', ml3)

    # tmp12: ml0, full, tmp68: ml2, full, L15: ml3, lower_triangular
    # tmp70 = (tmp12 L15^-T)
    trsm!('R', 'L', 'T', 'N', 1.0, ml3, ml0)

    # tmp68: ml2, full, L15: ml3, lower_triangular, tmp70: ml0, full
    ml4 = Array{Float64}(undef, 230)
    # tmp23 = (tmp70^T tmp68)
    gemv!('T', 1.0, ml0, ml2, 0.0, ml4)

    # L15: ml3, lower_triangular, tmp23: ml4, full
    # tmp24 = (L15^-T tmp23)
    trsv!('L', 'T', 'N', ml3, ml4)

    # tmp24: ml4, full
    # b = tmp24
    return (ml4)
end