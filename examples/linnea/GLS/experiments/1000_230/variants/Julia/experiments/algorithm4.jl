using LinearAlgebra.BLAS
using LinearAlgebra

function algorithm4(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    start::Float64 = 0.0
    finish::Float64 = 0.0
    Benchmarker.cachescrub()
    start = time_ns()

    # cost: 6.22e+08 FLOPs
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
    ml3 = Array{Float64}(undef, 230)
    # tmp21 = (tmp12^T tmp68)
    gemv!('T', 1.0, ml0, ml2, 0.0, ml3)

    # tmp12: ml0, full, tmp21: ml3, full
    ml4 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    syrk!('L', 'T', 1.0, ml0, 0.0, ml4)

    # tmp21: ml3, full, tmp14: ml4, symmetric_lower_triangular
    # (L15 L15^T) = tmp14
    LAPACK.potrf!('L', ml4)

    # tmp21: ml3, full, L15: ml4, lower_triangular
    # tmp23 = (L15^-1 tmp21)
    trsv!('L', 'N', 'N', ml4, ml3)

    # L15: ml4, lower_triangular, tmp23: ml3, full
    # tmp24 = (L15^-T tmp23)
    trsv!('L', 'T', 'N', ml4, ml3)

    # tmp24: ml3, full
    # b = tmp24

    finish = time_ns()
    return (tuple(ml3), (finish-start)*1e-9)
end