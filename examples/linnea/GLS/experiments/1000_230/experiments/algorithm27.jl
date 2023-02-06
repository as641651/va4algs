using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm27(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})

Compute
b = ((X^T M^-1 X)^-1 X^T M^-1 y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1000 x 230 with property FullRank.
- `ml1::Array{Float64,2}`: Matrix M of size 1000 x 1000 with property SPD.
- `ml2::Array{Float64,1}`: Vector y of size 1000.
"""                    
function algorithm27(ml0::Array{Float64,2}, ml1::Array{Float64,2}, ml2::Array{Float64,1})
    # cost: 6.5e+08 FLOPs
    # X: ml0, full, M: ml1, full, y: ml2, full
    # (L2 L2^T) = M
    stime0= time()
    LAPACK.potrf!('L', ml1)

    # X: ml0, full, y: ml2, full, L2: ml1, lower_triangular
    # tmp68 = (L2^-1 y)
    stime1= time()
    trsv!('L', 'N', 'N', ml1, ml2)

    # X: ml0, full, L2: ml1, lower_triangular, tmp68: ml2, full
    # tmp12 = (L2^-1 X)
    stime2= time()
    trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0)

    # tmp68: ml2, full, tmp12: ml0, full
    ml3 = Array{Float64}(undef, 230, 230)
    # tmp14 = (tmp12^T tmp12)
    stime3= time()
    syrk!('L', 'T', 1.0, ml0, 0.0, ml3)

    # tmp68: ml2, full, tmp12: ml0, full, tmp14: ml3, symmetric_lower_triangular
    for i = 1:230-1;
        view(ml3, i, i+1:230)[:] = view(ml3, i+1:230, i);
    end;
    # (Q16 R17) = tmp14
    stime4= time()
    ml3 = qr!(ml3)

    # tmp68: ml2, full, tmp12: ml0, full, Q16: ml3, QRfact_Q, R17: ml3, QRfact_R
    ml4 = Array{Float64}(undef, 230)
    # tmp21 = (tmp12^T tmp68)
    stime5= time()
    gemv!('T', 1.0, ml0, ml2, 0.0, ml4)

    # Q16: ml3, QRfact_Q, R17: ml3, QRfact_R, tmp21: ml4, full
    ml5 = Array(ml3.Q)
    ml6 = Array{Float64}(undef, 230)
    # tmp25 = (Q16^T tmp21)
    stime6= time()
    gemv!('T', 1.0, ml5, ml4, 0.0, ml6)

    # R17: ml3, QRfact_R, tmp25: ml6, full
    ml7 = ml3.R
    # tmp24 = (R17^-1 tmp25)
    stime7= time()
    trsv!('U', 'N', 'N', ml7, ml6)

    # tmp24: ml6, full
    # b = tmp24
    stime8= time()
    return (ml6) ,(stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7,stime8,)
end

function write_algorithm27_to_eventlog(io, id, stamps)
    write( io, string(id, ";", "LAPACK.potrf_3.33e+08;", "3.33e+08;", "(L2 L2^T) = M;", "LAPACK.potrf!('L', ml1);", string(stamps[1]), ";", string(stamps[2]), '
'  ))
    write( io, string(id, ";", "trsv_1e+06;", "1e+06;", "tmp68 = (L2^-1 y);", "trsv!('L', 'N', 'N', ml1, ml2);", string(stamps[2]), ";", string(stamps[3]), '
'  ))
    write( io, string(id, ";", "trsm_2.3e+08;", "2.3e+08;", "tmp12 = (L2^-1 X);", "trsm!('L', 'L', 'N', 'N', 1.0, ml1, ml0);", string(stamps[3]), ";", string(stamps[4]), '
'  ))
    write( io, string(id, ";", "syrk_5.29e+07;", "5.29e+07;", "tmp14 = (tmp12^T tmp12);", "syrk!('L', 'T', 1.0, ml0, 0.0, ml3);", string(stamps[4]), ";", string(stamps[5]), '
'  ))
    write( io, string(id, ";", "ml3 = qr_3.24e+07;", "3.24e+07;", "(Q16 R17) = tmp14;", "ml3 = qr!(ml3);", string(stamps[5]), ";", string(stamps[6]), '
'  ))
    write( io, string(id, ";", "gemv_4.6e+05;", "4.6e+05;", "tmp21 = (tmp12^T tmp68);", "gemv!('T', 1.0, ml0, ml2, 0.0, ml4);", string(stamps[6]), ";", string(stamps[7]), '
'  ))
    write( io, string(id, ";", "gemv_1.06e+05;", "1.06e+05;", "tmp25 = (Q16^T tmp21);", "gemv!('T', 1.0, ml5, ml4, 0.0, ml6);", string(stamps[7]), ";", string(stamps[8]), '
'  ))
    write( io, string(id, ";", "trsv_5.29e+04;", "5.29e+04;", "tmp24 = (R17^-1 tmp25);", "trsv!('U', 'N', 'N', ml7, ml6);", string(stamps[8]), ";", string(stamps[9]), '
'  ))
end