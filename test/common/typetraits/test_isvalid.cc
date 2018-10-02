#include <config.h>

#include <iostream>
#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {
namespace Test {

struct MyVector {
public:
    void resize(std::size_t i) {};

    void resize(std::size_t i) const {};

    void resize() {};
};

struct MyOtherVector {
    double resize;
};

const auto hasResize = isValid([](auto&& a) -> decltype(a.resize(std::size_t(1))) { });

// using the check function
template<class Vector, typename std::enable_if_t<decltype(hasResize.template check<Vector>())::value, int> = 0>
void resize(const Vector& v, std::size_t size)
{
    v.resize(size); std::cout << "-> resized resizeable vector! size: " << size << std::endl;
}

// using decltype + declval and operator ()
template<class Vector, typename std::enable_if_t<!decltype(hasResize(std::declval<Vector>()))::value, int> = 0>
void resize(const Vector& v, std::size_t size)
{
    std::cout << "-> Did not resize non-resizeable vector!" << std::endl;
}

// using trailing return type and operator ()
template<class Vector>
auto resize2(const Vector& v, std::size_t size)
-> typename std::enable_if_t<decltype(hasResize(v))::value, void>
{
    v.resize(size); std::cout << "-> resized resizeable vector! size: " << size << std::endl;
}

// using trailing return type and operator ()
template<class Vector>
auto resize2(const Vector& v, std::size_t size)
-> typename std::enable_if_t<!decltype(hasResize.template check<Vector>())::value, void>
{
    std::cout << "-> Did not resize non-resizeable vector!" << std::endl;
}

} // end namespace Test
} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux::Test;

    MyVector v0;
    MyOtherVector v1;

    static_assert(hasResize(v0), "Vector v0 doesn't have a resize function!");
    static_assert(!hasResize(v1), "False positive. v1 has no resize function!");

    resize(v0, 3);
    resize(v1, 3);
    resize2(v0, 10);
    resize2(v1, 10);

    return 0;
}
