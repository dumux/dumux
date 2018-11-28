#include <config.h>

#include <iostream>
#include <tuple>

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
}

} // end namespace Dumux

int main(int argc, char* argv[]) try
{
    using namespace Dumux;

    runTest<Dune::MultiTypeBlockVector>();
    runTest<std::tuple>();

    return 0;

}
catch (const Dune::Exception& e)
{
    std::cout << e << std::endl;
    return 1;
}
catch (...)
{
    std::cout << "Unknown exception thrown!" << std::endl;
    return 1;
}
