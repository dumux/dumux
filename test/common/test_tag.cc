#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/common/tag.hh>

namespace Dumux {

struct T1 : public Utility::Tag<T1> {};
struct T2 : public Utility::Tag<T2> { std::string name() const { return "customname_t2"; }};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    T1 t11, t12;
    T2 t2;

    static_assert(t11 == t12, "Same tags not equal!");
    static_assert(t11 != t2, "Different tags should not be equal!");

    {
        std::stringstream s; s << t11;
        if (s.str() != "T1")
            DUNE_THROW(Dune::Exception, "Wrong name: " << s.str() << ", expected: T1.");
    }{
        std::stringstream s; s << t2;
        if (s.str() != "customname_t2")
            DUNE_THROW(Dune::Exception, "Wrong name: " << s.str() << ", expected: customname_t2.");
    }

    return 0;
}
