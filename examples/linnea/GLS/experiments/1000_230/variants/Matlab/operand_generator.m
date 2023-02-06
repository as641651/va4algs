function [out] = operand_generator()
    import MatrixGenerator.*;
    out{ 1 } = generate([1000,230], Shape.General(), Properties.Random([-1, 1]));
    out{ 2 } = generate([1000,1000], Shape.Symmetric(), Properties.SPD());
    out{ 3 } = generate([1000,1], Shape.General(), Properties.Random([-1, 1]));
end