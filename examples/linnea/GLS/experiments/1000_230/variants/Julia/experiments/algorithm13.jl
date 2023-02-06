using LinearAlgebra.BLAS
using LinearAlgebra

function algorithm13(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    start::Float64 = 0.0
    finish::Float64 = 0.0
    Benchmarker.cachescrub()
    start = time_ns()

    # cost: 6.23e+08 FLOPs
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
    # tmp68 = (L2^-1 y)
    trsv!('L', 'N', 'N', ml1, ml2)

    # tmp56: ml3, full, L2: ml1, lower_triangular, tmp12: ml0, full, tmp68: ml2, full
    ml4 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    syrk!('L', 'T', 1.0, ml0, 0.0, ml4)

    # tmp56: ml3, full, L2: ml1, lower_triangular, tmp68: ml2, full, tmp14: ml4, symmetric_lower_triangular
    # tmp71 = (L2^-T tmp68)
    trsv!('L', 'T', 'N', ml1, ml2)

    # tmp56: ml3, full, tmp14: ml4, symmetric_lower_triangular, tmp71: ml2, full
    # (L15 L15^T) = tmp14
    LAPACK.potrf!('L', ml4)

    # tmp56: ml3, full, tmp71: ml2, full, L15: ml4, lower_triangular
    ml5 = Array{Float64}(undef, 230)
    # tmp21 = (tmp56 tmp71)
    gemv!('N', 1.0, ml3, ml2, 0.0, ml5)

    # L15: ml4, lower_triangular, tmp21: ml5, full
    # tmp23 = (L15^-1 tmp21)
    trsv!('L', 'N', 'N', ml4, ml5)

    # L15: ml4, lower_triangular, tmp23: ml5, full
    # tmp24 = (L15^-T tmp23)
    trsv!('L', 'T', 'N', ml4, ml5)

    # tmp24: ml5, full
    # b = tmp24

    finish = time_ns()
    return (tuple(ml5), (finish-start)*1e-9)
end