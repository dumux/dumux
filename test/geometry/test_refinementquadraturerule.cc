#include <config.h>

#include <iostream>
#include <array>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dumux/geometry/refinementquadraturerule.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // make a simplex
    const auto simplex = Dune::AffineGeometry<double, 2, 3>(
        Dune::GeometryTypes::simplex(2),
        std::array<Dune::FieldVector<double, 3>, 3>{{{0.0, 0.0, 1.0}, {2.0, 3.0, 4.0}, {3.0, -1.0, 1.0}}}
    );

    // make quadrature rule
    const auto quad = RefinementQuadratureRule<double, 2>(Dune::GeometryTypes::simplex(2), /*refinementlevels=*/3);

    // integrate volume
    double volume = 0.0;
    for (const auto& qp : quad)
        volume += qp.weight()*simplex.integrationElement(qp.position());

    if (Dune::FloatCmp::ne(volume, simplex.volume()))
        DUNE_THROW(Dune::Exception, "Integration yields wrong simplex volume: " << volume << ", should be " << simplex.volume());

    return 0;
}
