using LinearAlgebra.BLAS
using LinearAlgebra

function algorithm55(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    start::Float64 = 0.0
    finish::Float64 = 0.0
    Benchmarker.cachescrub()
    start = time_ns()

    # cost: 6.59e+08 FLOPs
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
    ml4 = Array{Float64}(undef, 230)
    # (Z18 W19 Z18^T) = tmp14
    ml4, ml3 = LAPACK.syev!('V', 'L', ml3)

    # tmp12: ml0, full, tmp68: ml2, full, W19: ml4, diagonal_vector, Z18: ml3, full
    ml5 = Array{Float64}(undef, 230)
    # tmp21 = (tmp12^T tmp68)
    gemv!('T', 1.0, ml0, ml2, 0.0, ml5)

    # W19: ml4, diagonal_vector, Z18: ml3, full, tmp21: ml5, full
    ml6 = Array{Float64}(undef, 230)
    # tmp268 = (Z18^T tmp21)
    gemv!('T', 1.0, ml3, ml5, 0.0, ml6)

    # W19: ml4, diagonal_vector, Z18: ml3, full, tmp268: ml6, full
    # tmp269 = (W19^-1 tmp268)
    ml6 ./= ml4

    # Z18: ml3, full, tmp269: ml6, full
    ml7 = Array{Float64}(undef, 230)
    # tmp24 = (Z18 tmp269)
    gemv!('N', 1.0, ml3, ml6, 0.0, ml7)

    # tmp24: ml7, full
    # b = tmp24

    finish = time_ns()
    return (tuple(ml7), (finish-start)*1e-9)
end