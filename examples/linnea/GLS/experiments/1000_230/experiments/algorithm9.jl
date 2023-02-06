using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm9(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm9(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.23e+08 FLOPs
    # X: ml0, full, M: ml1, full, y: ml2, full
    ml3 = Array{Float64}(undef, 230, 1000)
    # tmp56 = X^T
    stime0= time()
    transpose!(ml3, ml0)

    # X: ml0, full, M: ml1, full, y: ml2, full, tmp56: ml3, full
    # (L2 L2^T) = M
    stime1= time()
    LAPACK.potrf!('L', ml1)

    # X: ml0, full, y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular
    # tmp12 = (L2^-1 X)
    stime2= time()
    trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0)

    # y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular, tmp12: ml0, full
    ml4 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    stime3= time()
    syrk!('L', 'T', 1.0, ml0, 0.0, ml4)

    # y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular, tmp14: ml4, symmetric_lower_triangular
    # (L15 L15^T) = tmp14
    stime4= time()
    LAPACK.potrf!('L', ml4)

    # y: ml2, full, tmp56: ml3, full, L2: ml1, lower_triangular, L15: ml4, lower_triangular
    # tmp68 = (L2^-1 y)
    stime5= time()
    trsv!('L', 'N', 'N', ml1, ml2)

    # tmp56: ml3, full, L2: ml1, lower_triangular, L15: ml4, lower_triangular, tmp68: ml2, full
    # tmp71 = (L2^-T tmp68)
    stime6= time()
    trsv!('L', 'T', 'N', ml1, ml2)

    # tmp56: ml3, full, L15: ml4, lower_triangular, tmp71: ml2, full
    ml5 = Array{Float64}(undef, 230)
    # tmp21 = (tmp56 tmp71)
    stime7= time()
    gemv!('N', 1.0, ml3, ml2, 0.0, ml5)

    # L15: ml4, lower_triangular, tmp21: ml5, full
    # tmp23 = (L15^-1 tmp21)
    stime8= time()
    trsv!('L', 'N', 'N', ml4, ml5)

    # L15: ml4, lower_triangular, tmp23: ml5, full
    # tmp24 = (L15^-T tmp23)
    stime9= time()
    trsv!('L', 'T', 'N', ml4, ml5)

    # tmp24: ml5, full
    # b = tmp24
    stime10= time()
    return (ml5) ,(stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7,stime8,stime9,stime10,)
end

function write_algorithm9_to_eventlog(io, id, stamps)
    write( io, string(id, ";", "transpose_1;", "1;", "tmp56 = X^T;", "transpose!(ml3, ml0);", string(stamps[1]), ";", string(stamps[2]), '
'  ))
    write( io, string(id, ";", "LAPACK.potrf_3.33e+08;", "3.33e+08;", "(L2 L2^T) = M;", "LAPACK.potrf!('L', ml1);", string(stamps[2]), ";", string(stamps[3]), '
'  ))
    write( io, string(id, ";", "trsm_2.3e+08;", "2.3e+08;", "tmp12 = (L2^-1 X);", "trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0);", string(stamps[3]), ";", string(stamps[4]), '
'  ))
    write( io, string(id, ";", "syrk_5.29e+07;", "5.29e+07;", "tmp14 = (tmp12^T tmp12);", "syrk!('L', 'T', 1.0, ml0, 0.0, ml4);", string(stamps[4]), ";", string(stamps[5]), '
'  ))
    write( io, string(id, ";", "LAPACK.potrf_4.06e+06;", "4.06e+06;", "(L15 L15^T) = tmp14;", "LAPACK.potrf!('L', ml4);", string(stamps[5]), ";", string(stamps[6]), '
'  ))
    write( io, string(id, ";", "trsv_1e+06;", "1e+06;", "tmp68 = (L2^-1 y);", "trsv!('L', 'N', 'N', ml1, ml2);", string(stamps[6]), ";", string(stamps[7]), '
'  ))
    write( io, string(id, ";", "trsv_1e+06;", "1e+06;", "tmp71 = (L2^-T tmp68);", "trsv!('L', 'T', 'N', ml1, ml2);", string(stamps[7]), ";", string(stamps[8]), '
'  ))
    write( io, string(id, ";", "gemv_4.6e+05;", "4.6e+05;", "tmp21 = (tmp56 tmp71);", "gemv!('N', 1.0, ml3, ml2, 0.0, ml5);", string(stamps[8]), ";", string(stamps[9]), '
'  ))
    write( io, string(id, ";", "trsv_5.29e+04;", "5.29e+04;", "tmp23 = (L15^-1 tmp21);", "trsv!('L', 'N', 'N', ml4, ml5);", string(stamps[9]), ";", string(stamps[10]), '
'  ))
    write( io, string(id, ";", "trsv_5.29e+04;", "5.29e+04;", "tmp24 = (L15^-T tmp23);", "trsv!('L', 'T', 'N', ml4, ml5);", string(stamps[10]), ";", string(stamps[11]), '
'  ))
end