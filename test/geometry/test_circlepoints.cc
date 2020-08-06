#include <config.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>

#include <dumux/multidomain/embedded/circlepoints.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;

    using Point = Dune::FieldVector<double, 3>;
    const Point center(0.0);
    const Point normal({1.0, 1.0, 0.0}); // doesn't need to be a unit normal vector
    const double radius = 3.0;
    static constexpr std::size_t numPoints = 100;

    // make a circle and and ellipse with equal major radii -> should yield the same points
    const auto circle = EmbeddedCoupling::circlePoints(center, normal, radius, numPoints);
    const auto checkRadius = [&](const auto& p){ return Dune::FloatCmp::eq(p.two_norm(), radius, 1e-7); };
    if (!std::all_of(circle.begin(), circle.end(), checkRadius))
        DUNE_THROW(Dune::Exception, "Not all point have the same radius!");

    Point unitNormal = normal;
    unitNormal /= unitNormal.two_norm();
    bool correctOrientation = true;
    for (int i = 1; i < circle.size(); ++i)
    {
        auto n = crossProduct(circle[i]-center, circle[i-1]-center);
        n /= n.two_norm();
        correctOrientation &= Dune::FloatCmp::eq(std::abs(n*unitNormal), 1.0, 1e-10);
    }
    if (!correctOrientation)
        DUNE_THROW(Dune::Exception, "Not all points have the correct orientation!");

    if (circle.size() != numPoints)
        DUNE_THROW(Dune::Exception, "Wrong number of points!");

    return 0;
}
