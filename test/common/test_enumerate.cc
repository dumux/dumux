#include <config.h>

#include <vector>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dumux/common/enumerate.hh>

int main(int argc, char* argv[]) {

    using namespace Dumux;

    std::vector<std::size_t> vec(20);
    std::iota(vec.begin(), vec.end(), 0);

    for (const auto& [i, item] : enumerate(vec))
        if (i != item)
            DUNE_THROW(Dune::Exception, "Wrong index!");

    return 0;
}
