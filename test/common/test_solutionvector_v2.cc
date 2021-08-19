#include <utility>
#include <iostream>
#include <algorithm>

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/vectoradapter.hh>
#include <dumux/common/vectorwithstate.hh>

#ifndef VECTORSIZE
#define VECTORSIZE 2
#endif

// functions to make instances of the vectors used in the tests
namespace TestFactory {

namespace Detail {

struct HasResize
{
    template<class V>
    auto operator()(const V& v) -> decltype(std::declval<V>().resize(0)) {}
};

template<class Vector>
inline constexpr bool hasResize = decltype(
    Dumux::isValid(HasResize())(std::declval<Vector>())
)::value;

} // end namespace Detail

template<typename PriVars, typename Value>
void fillPriVars(PriVars& p, Value v)
{ std::fill(p.begin(), p.end(), v); }

template<typename Vector>
Vector makeTestVector()
{
    Vector v;
    if constexpr (Detail::hasResize<Vector>)
        v.resize(VECTORSIZE);
    for (auto& p : v)
        fillPriVars(p, 1);
    return v;
}

template<typename MultiTypeVector>
MultiTypeVector makeMultiTypeTestVector()
{
    MultiTypeVector v;

    using namespace Dune::Hybrid;
    static constexpr std::size_t N = MultiTypeVector::size();
    forEach(std::make_index_sequence<N>{}, [&] (auto i) {
        using V = std::decay_t<decltype(v[Dune::index_constant<i>()])>;
        v[Dune::index_constant<i>()] = makeTestVector<V>();
    });

    return v;
}

template<typename Vector>
Vector makeTestVectorWithState()
{
    auto v = makeTestVector<Vector>();
    for (auto& p : v)
        p.setState(1);
    return v;
}

template<typename MultiTypeVector>
MultiTypeVector makeMultiTypeTestVectorWithState()
{
    MultiTypeVector v;

    using namespace Dune::Hybrid;
    static constexpr std::size_t N = MultiTypeVector::size();
    forEach(std::make_index_sequence<N>{}, [&] (auto i) {
        using V = std::decay_t<decltype(v[Dune::index_constant<i>()])>;
        v[Dune::index_constant<i>()] = makeTestVectorWithState<V>();
    });

    return v;
}

} // end namespace TestFactory


// functions to test equality of privars & states in the test vectors
namespace TestEQCheck {

template<typename V1, typename V2>
bool priVars(const V1& v1, const V2& v2)
{
    const auto equalPriVars = [] (const auto& v1, const auto& v2) {
        using std::abs;
        return std::equal(v1.begin(), v1.end(), v2.begin(), v2.end(),
                          [] (const auto& a, const auto& b) {
                              return abs(a-b) < 1e-6;
                          });
    };
    return std::equal(v1.begin(), v1.end(),
                      v2.begin(), v2.end(),
                      equalPriVars);
}

template<typename V1, typename V2>
bool states(const V1& v1, const V2& v2)
{
    const auto equalState = [] (const auto& v1, const auto& v2) {
        return v1.state() == v2.state();
    };
    return std::equal(v1.begin(), v1.end(),
                      v2.begin(), v2.end(),
                      equalState);
}

template<typename V1, typename V2>
bool statesAndPriVars(const V1& v1, const V2& v2)
{ return states(v1, v2) && priVars(v1, v2); }

template<typename V1, typename V2>
bool multiTypePriVars(const V1& v1, const V2& v2)
{
    if (v1.size() != v2.size())
        return false;

    bool eq = true;
    using namespace Dune::Hybrid;
    static constexpr std::size_t N = V1::size();
    forEach(std::make_index_sequence<N>{}, [&] (auto i) {
        eq = eq && priVars(v1[Dune::index_constant<i>{}],
                           v2[Dune::index_constant<i>{}]);
    });
    return eq;
}

template<typename V1, typename V2>
bool multiTypeStates(const V1& v1, const V2& v2)
{
    if (v1.size() != v2.size())
        return false;

    bool eq = true;
    using namespace Dune::Hybrid;
    static constexpr std::size_t N = V1::size();
    forEach(std::make_index_sequence<N>{}, [&] (auto i) {
        eq = eq && states(v1[Dune::index_constant<i>{}],
                          v2[Dune::index_constant<i>{}]);
    });
    return eq;
}

template<typename V1, typename V2>
bool multiTypeStatesAndPriVars(const V1& v1, const V2& v2)
{ return multiTypeStates(v1, v2) && multiTypePriVars(v1, v2); }

} // end namespace TestEQCheck


// the actual test functions
namespace Test {

template<typename Vector>
void printTestMessage()
{
    std::cout << "Testing with \"" << Dune::className<Vector>() << "\"" << std::endl;
}

template<typename Vector, typename EqualityCheck>
void operators(const Vector& v, const EqualityCheck& equalityCheck)
{
    auto copy = v;
    if (!equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "Copy failed");

    copy *= 2;
    if (equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "multiply failed");

    copy /= 2;
    if (!equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "division failed");

    copy += copy;
    if (equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "addition failed");

    copy -= copy;
    if (equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "subtraction failed");
}

template<typename Vector, typename EqualityCheck>
void states(const Vector& v, const EqualityCheck& equalityCheck)
{
    auto copy = v;
    if (!equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "Copy failed");

    const auto store = copy[0].state();
    copy[0].setState(store + store);
    if (equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "State modification failed");

    copy[0].setState(store);
    if (!equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "State reset failed");
}

template<typename Vector, typename EqualityCheck>
void multiTypeStates(const Vector& v, const EqualityCheck& equalityCheck)
{
    auto copy = v;
    if (!equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "Copy failed");

    using namespace Dune::Indices;
    const auto store = copy[_0][0].state();
    copy[_0][0].setState(store + store);
    if (equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "State modification failed");

    copy[_0][0].setState(store);
    if (!equalityCheck(copy, v))
        DUNE_THROW(Dune::InvalidStateException, "State reset failed");
}

template<typename Vector>
void vector()
{
    printTestMessage<Vector>();
    const auto v = TestFactory::makeTestVector<Vector>();
    const auto eq = [] (const auto& a, const auto& b) {
        return TestEQCheck::priVars(a, b);
    };
    operators(v, eq);
}

template<typename Vector>
void vectorWithState()
{
    printTestMessage<Vector>();
    const auto v = TestFactory::makeTestVectorWithState<Vector>();
    const auto eq = [] (const auto& a, const auto& b) {
        return TestEQCheck::statesAndPriVars(a, b);
    };
    operators(v, eq);
    states(v, eq);
}

template<typename Vector>
void multiTypeVector()
{
    printTestMessage<Vector>();
    const auto v = TestFactory::makeMultiTypeTestVector<Vector>();
    const auto eq = [] (const auto& a, const auto& b) {
        return TestEQCheck::multiTypePriVars(a, b);
    };
    operators(v, eq);
}

template<typename Vector>
void multiTypeVectorWithState()
{
    printTestMessage<Vector>();
    const auto v = TestFactory::makeMultiTypeTestVectorWithState<Vector>();
    const auto eq = [] (const auto& a, const auto& b) {
        return TestEQCheck::multiTypeStatesAndPriVars(a, b);
    };
    operators(v, eq);
    multiTypeStates(v, eq);
}

} // end namespace Test

int main()
{
    using Scalar = int;
    using State = int;
    using RawPriVars = Dune::FieldVector<Scalar, VECTORSIZE>;
    using PriVarsWithState = Dumux::SwitchablePrimaryVariables<RawPriVars, State>;

    {
        using V = Dune::BlockVector<RawPriVars>;
        using MTV = Dune::MultiTypeBlockVector<V, V>;
        Test::vector<V>();
        Test::multiTypeVector<MTV>();
    }
    {
        using V = Dune::BlockVector<PriVarsWithState>;
        using MTV = Dune::MultiTypeBlockVector<V, V>;
        Test::vectorWithState<V>();
        Test::multiTypeVectorWithState<MTV>();
    }
    {
        using BV = Dune::BlockVector<RawPriVars>;
        using MTV = Dune::MultiTypeBlockVector<BV, BV>;
        Test::vector<Dumux::VectorAdapter<BV>>();
        Test::multiTypeVector<Dumux::VectorAdapter<MTV>>();

        // test adapter with ref & const ref
        using FV = Dune::FieldVector<RawPriVars, 2>;
        FV v; v[0] = 1; v[1] = 2;
        Dumux::VectorViewAdapter<FV> constWrap(v);
        Dumux::VectorProxyAdapter<FV> wrap(v);

        if (wrap[0][0] != v[0][0] || wrap[0][1] != v[0][1])
            DUNE_THROW(Dune::InvalidStateException, "Wrong wrapped value");
        if (wrap[1][0] != v[1][0] || wrap[1][1] != v[1][1])
            DUNE_THROW(Dune::InvalidStateException, "Wrong wrapped value");
        if (constWrap[0][0] != v[0][0] || constWrap[0][1] != v[0][1])
            DUNE_THROW(Dune::InvalidStateException, "Wrong const wrapped value");
        if (constWrap[1][0] != v[1][0] || constWrap[1][1] != v[1][1])
            DUNE_THROW(Dune::InvalidStateException, "Wrong const wrapped value");

        // value modification should affect both vectors
        const auto tmp = wrap[0][0]*2;
        wrap[0][0] = tmp;
        if (v[0][0] != wrap[0][0])
            DUNE_THROW(Dune::InvalidStateException, "Value modification failed");
    }
    {
        using BV = Dune::BlockVector<RawPriVars>;
        using States = std::vector<int>;
        Dumux::VectorWithState<BV, States> v;

        v.resize(2);
        v[0] = 1; v[1] = 2;
        if (v[0][0] != 1 && v[0][1] != 2)
            DUNE_THROW(Dune::InvalidStateException, "Value set failed");

        v[0].setState(10); v[1].setState(12);
        if (v[0].state() != 10 && v[1].state() != 12)
            DUNE_THROW(Dune::InvalidStateException, "State set failed");
    }

    return 0;
}
