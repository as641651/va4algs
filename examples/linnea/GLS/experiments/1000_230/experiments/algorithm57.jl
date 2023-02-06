using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm57(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm57(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.59e+08 FLOPs
    # X: ml0, full, M: ml1, full, y: ml2, full
    # (L2 L2^T) = M
    stime0= time()
    LAPACK.potrf!('L', ml1)

    # X: ml0, full, y: ml2, full, L2: ml1, lower_triangular
    # tmp12 = (L2^-1 X)
    stime1= time()
    trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0)

    # y: ml2, full, L2: ml1, lower_triangular, tmp12: ml0, full
    # tmp68 = (L2^-1 y)
    stime2= time()
    trsv!('L', 'N', 'N', ml1, ml2)

    # tmp12: ml0, full, tmp68: ml2, full
    ml3 = Array{Float64}(undef, 230)
    # tmp21 = (tmp12^T tmp68)
    stime3= time()
    gemv!('T', 1.0, ml0, ml2, 0.0, ml3)

    # tmp12: ml0, full, tmp21: ml3, full
    ml4 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    stime4= time()
    syrk!('L', 'T', 1.0, ml0, 0.0, ml4)

    # tmp21: ml3, full, tmp14: ml4, symmetric_lower_triangular
    ml5 = Array{Float64}(undef, 230)
    # (Z18 W19 Z18^T) = tmp14
    stime5= time()
    ml5, ml4 = LAPACK.syev!('V', 'L', ml4)

    # tmp21: ml3, full, W19: ml5, diagonal_vector, Z18: ml4, full
    ml6 = Array{Float64}(undef, 230)
    # tmp268 = (Z18^T tmp21)
    stime6= time()
    gemv!('T', 1.0, ml4, ml3, 0.0, ml6)

    # W19: ml5, diagonal_vector, Z18: ml4, full, tmp268: ml6, full
    # tmp269 = (W19^-1 tmp268)
    ml6 ./= ml5

    # Z18: ml4, full, tmp269: ml6, full
    ml7 = Array{Float64}(undef, 230)
    # tmp24 = (Z18 tmp269)
    stime7= time()
    gemv!('N', 1.0, ml4, ml6, 0.0, ml7)

    # tmp24: ml7, full
    # b = tmp24
    stime8= time()
    return (ml7) ,(stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7,stime8,)
end

function write_algorithm57_to_eventlog(io, id, stamps)
    write( io, string(id, ";", "LAPACK.potrf_3.33e+08;", "3.33e+08;", "(L2 L2^T) = M;", "LAPACK.potrf!('L', ml1);", string(stamps[1]), ";", string(stamps[2]), '
'  ))
    write( io, string(id, ";", "trsm_2.3e+08;", "2.3e+08;", "tmp12 = (L2^-1 X);", "trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0);", string(stamps[2]), ";", string(stamps[3]), '
'  ))
    write( io, string(id, ";", "trsv_1e+06;", "1e+06;", "tmp68 = (L2^-1 y);", "trsv!('L', 'N', 'N', ml1, ml2);", string(stamps[3]), ";", string(stamps[4]), '
'  ))
    write( io, string(id, ";", "gemv_4.6e+05;", "4.6e+05;", "tmp21 = (tmp12^T tmp68);", "gemv!('T', 1.0, ml0, ml2, 0.0, ml3);", string(stamps[4]), ";", string(stamps[5]), '
'  ))
    write( io, string(id, ";", "syrk_5.29e+07;", "5.29e+07;", "tmp14 = (tmp12^T tmp12);", "syrk!('L', 'T', 1.0, ml0, 0.0, ml4);", string(stamps[5]), ";", string(stamps[6]), '
'  ))
    write( io, string(id, ";", "ml5, ml4 = LAPACK.syev_4.06e+07;", "4.06e+07;", "(Z18 W19 Z18^T) = tmp14;", "ml5, ml4 = LAPACK.syev!('V', 'L', ml4);", string(stamps[6]), ";", string(stamps[7]), '
'  ))
    write( io, string(id, ";", "gemv_1.06e+05;", "1.06e+05;", "tmp268 = (Z18^T tmp21);", "gemv!('T', 1.0, ml4, ml3, 0.0, ml6);", string(stamps[7]), ";", string(stamps[8]), '
'  ))
    write( io, string(id, ";", "gemv_230;", "230;", "tmp269 = (W19^-1 tmp268);", "gemv!('N', 1.0, ml4, ml6, 0.0, ml7);", string(stamps[8]), ";", string(stamps[9]), '
'  ))
end