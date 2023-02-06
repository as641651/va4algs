using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm39(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm39(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.51e+08 FLOPs
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
    for i = 1:230-1;
        view(ml4, i, i+1:230)[:] = view(ml4, i+1:230, i);
    end;
    # (Q16 R17) = tmp14
    ml4 = qr!(ml4)

    # tmp56: ml3, full, L2: ml1, lower_triangular, tmp68: ml2, full, Q16: ml4, QRfact_Q, R17: ml4, QRfact_R
    # tmp71 = (L2^-T tmp68)
    trsv!('L', 'T', 'N', ml1, ml2)

    # tmp56: ml3, full, Q16: ml4, QRfact_Q, R17: ml4, QRfact_R, tmp71: ml2, full
    ml5 = Array{Float64}(undef, 230)
    # tmp21 = (tmp56 tmp71)
    gemv!('N', 1.0, ml3, ml2, 0.0, ml5)

    # Q16: ml4, QRfact_Q, R17: ml4, QRfact_R, tmp21: ml5, full
    ml6 = Array(ml4.Q)
    ml7 = Array{Float64}(undef, 230)
    # tmp25 = (Q16^T tmp21)
    gemv!('T', 1.0, ml6, ml5, 0.0, ml7)

    # R17: ml4, QRfact_R, tmp25: ml7, full
    ml8 = ml4.R
    # tmp24 = (R17^-1 tmp25)
    trsv!('U', 'N', 'N', ml8, ml7)

    # tmp24: ml7, full
    # b = tmp24
    return (ml7)
end