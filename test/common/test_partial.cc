#include <config.h>

#include <iostream>
#include <tuple>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/partial.hh>

namespace Dumux {

template<template<class...> class T>
void runTest()
{
    using Block1 = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using Block2 = Dune::BlockVector<Dune::FieldVector<double, 3>>;

    Block1 a, b;
    Block2 c;

    a = {1, 3, 4};
    b = {4, 5, 6, 7};
    c = {{1, 2, 3}, {1, 2, 3}};

    using namespace Dune::Indices;

    T<Block1, Block1, Block2> m;

    std::cout << "Testing " << Dune::className(m) << '\n' << std::endl;

    std::get<0>(m) = a;
    std::get<1>(m) = b;
    std::get<2>(m) = c;

    auto p = partial(m, _0, _2);
    p = partial(m, std::make_tuple(_0, _2));

    if (!std::is_same<T<Block1&, Block2&>, std::decay_t<decltype(p)>>::value)
        DUNE_THROW(Dune::Exception, "Dumux::partial() returned wrong type: " << Dune::className(p));

    std::get<1>(p)[0][0] = 5.0;

    if (!Dune::FloatCmp::eq(std::get<2>(m)[0][0], 5.0))
        DUNE_THROW(Dune::Exception, "Modifying referenced partial vector failed! (m = " << std::get<2>(m)[0][0] << ", p = " << std::get<1>(p)[0][0] << ")");

    const auto mConst = m;
    auto pConst = partial(mConst, _0, _2);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::get<0>(pConst))>>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::get<1>(pConst))>>);

    if (!Dune::FloatCmp::eq(std::get<0>(mConst)[0][0], std::get<0>(pConst)[0][0]))
        DUNE_THROW(Dune::Exception, "Values differ! (mConst = " << std::get<0>(mConst)[0][0] << ", pConst = " << std::get<0>(pConst)[0][0] << ")");

    const auto& mConstRef = m;
    auto pConstRef = partial(mConstRef, _0, _2);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::get<0>(pConstRef))>>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::get<1>(pConstRef))>>);

    std::get<2>(m)[0][0] = 10.0;

    if (!Dune::FloatCmp::eq(std::get<1>(pConstRef)[0][0], 10.0))
        DUNE_THROW(Dune::Exception, "Modifying referenced partial vector failed! (m = " << std::get<2>(m)[0][0] << ", pConstRef = " << std::get<1>(pConstRef)[0][0] << ")");
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    runTest<Dune::MultiTypeBlockVector>();
    runTest<std::tuple>();

    return 0;
}
