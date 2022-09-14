#include <dune/common/exceptions.hh>
#include <dumux/experimental/new_assembly/dumux/common/size.hh>

int main()
{
    using Dumux::Size::dynamic;
    using Dumux::Size::isRepresentableBy;
    if (dynamic*dynamic != dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator*");
    if (dynamic+1 != dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator+");
    if (dynamic*1 != dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator*int{}");
    if (Dumux::Size::pow<dynamic, 2>() != dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected power result");

    if (dynamic != dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Expected dynamic sizes to always compare equal");
    if (!(dynamic == dynamic))
        DUNE_THROW(Dune::InvalidStateException, "Expected dynamic sizes to always compare equal");
    if (dynamic < dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator<");
    if (dynamic > dynamic)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator>");

    if (!isRepresentableBy(dynamic, dynamic))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected isRepresentable(dynamic, dynamic) result");
    if (isRepresentableBy(int{}, dynamic))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected isRepresentable(int, dynamic) result");
    if (!isRepresentableBy(dynamic, int{}))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected isRepresentable(dynamic, int) result");
    if (!isRepresentableBy(int{3}, int{2}))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected isRepresentable(3, 2) result");
    if (!isRepresentableBy(int{3}, int{3}))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected isRepresentable(3, 3) result");
    if (isRepresentableBy(int{2}, int{3}))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected isRepresentable(2, 3) result");

    return 0;
}
