struct naive_eigen
{
template<typename Type_X, typename Type_M, typename Type_y>
decltype(auto) operator()(Type_X && X, Type_M && M, Type_y && y)
{
    auto b = (((X).transpose()*(M).inverse()*X).inverse()*(X).transpose()*(M).inverse()*y).eval();

    typedef std::remove_reference_t<decltype(b)> return_t;
    return return_t(b);                         
}
};