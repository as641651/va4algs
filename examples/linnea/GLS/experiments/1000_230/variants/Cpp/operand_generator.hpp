#include <generator/generator.hpp>

template<typename Gen>
decltype(auto) operand_generator(Gen && gen)
{
    auto X = gen.generate({1000,230}, generator::property::random{}, generator::shape::not_square{});
    auto M = gen.generate({1000,1000}, generator::shape::self_adjoint{}, generator::property::spd{});
    auto y = gen.generate({1000,1}, generator::property::random{}, generator::shape::col_vector{}, generator::shape::not_square{});
    return std::make_tuple(X, M, y);
}