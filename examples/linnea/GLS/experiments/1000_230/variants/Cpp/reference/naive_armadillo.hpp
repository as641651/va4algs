struct naive_armadillo
{
template<typename Type_X, typename Type_M, typename Type_y>
decltype(auto) operator()(Type_X && X, Type_M && M, Type_y && y)
{
    auto b = (((X).t()*arma::inv_sympd(M)*X).i()*(X).t()*arma::inv_sympd(M)*y).eval();

    typedef std::remove_reference_t<decltype(b)> return_t;
    return return_t(b);                         
}
};