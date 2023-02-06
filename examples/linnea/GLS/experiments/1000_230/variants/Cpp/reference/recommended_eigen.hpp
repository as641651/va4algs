struct recommended_eigen
{
template<typename Type_X, typename Type_M, typename Type_y>
decltype(auto) operator()(Type_X && X, Type_M && M, Type_y && y)
{
    auto b = (( ((X).transpose()*( (M).llt().solve(X) )).partialPivLu().solve((X).transpose()) )*( (M).llt().solve(y) )).eval();

    typedef std::remove_reference_t<decltype(b)> return_t;
    return return_t(b);                         
}
};