using LinearAlgebra.BLAS
using LinearAlgebra

function recommended(X::Array{Float64,2}, M::Symmetric{Float64,Array{Float64,2}}, y::Array{Float64,1})
    start::Float64 = 0.0
    finish::Float64 = 0.0
    Benchmarker.cachescrub()
    start = time_ns()

    b = ((transpose(X)*((M)\X))\transpose(X))*((M)\y);

    finish = time_ns()
    return (tuple(b), (finish-start)*1e-9)
end