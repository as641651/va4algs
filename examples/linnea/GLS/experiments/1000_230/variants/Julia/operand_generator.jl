using MatrixGenerator

function operand_generator()
    X::Array{Float64,2} = generate((1000,230), [Shape.General, Properties.Random(-1, 1)])
    M::Symmetric{Float64,Array{Float64,2}} = generate((1000,1000), [Shape.Symmetric, Properties.SPD])
    y::Array{Float64,1} = generate((1000,1), [Shape.General, Properties.Random(-1, 1)])
    return (X, M, y,)
end