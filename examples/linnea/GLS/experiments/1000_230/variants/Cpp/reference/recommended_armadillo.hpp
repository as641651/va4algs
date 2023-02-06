struct recommended_armadillo
{
template<typename Type_X, typename Type_M, typename Type_y>
decltype(auto) operator()(Type_X && X, Type_M && M, Type_y && y)
{
    auto b = (arma::solve((X).t()*arma::solve(M, X, arma::solve_opts::fast), (X).t(), arma::solve_opts::fast)*arma::solve(M, y, arma::solve_opts::fast)).eval();

    typedef std::remove_reference_t<decltype(b)> return_t;
    return return_t(b);                         
}
};