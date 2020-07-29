#include <config.h>

#include <iostream>
#include <iomanip>
#include <cmath>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/multidomain/embedded/cylinderintegration.hh>
#include "transformation.hh"

namespace Dumux {

template<class Scalar>
void checkVolumeIntegral(const Scalar volIntegralComputed, const Scalar volIntegralExact)
{
    std::cout << std::setprecision(7) << std::scientific
              << "Volume integral: " << volIntegralComputed << " (computed), "
                                     << volIntegralExact << " (exact)" << "\n";
    if (Dune::FloatCmp::ne(volIntegralComputed, volIntegralExact, 1e-7))
        DUNE_THROW(Dune::Exception, "Volume integral is not exact!");
}

template<class Scalar>
void checkAreaIntegral(const Scalar areaIntegralComputed, const Scalar areaIntegralExact)
{
    std::cout << std::setprecision(7) << std::scientific
              << "Area integral: " << areaIntegralComputed << " (computed), "
                                   << areaIntegralExact << " (exact)" << "\n";
    if (Dune::FloatCmp::ne(areaIntegralComputed, areaIntegralExact, 1e-7))
        DUNE_THROW(Dune::Exception, "Area integral is not exact!");
}

template<class GlobalPosition>
void checkCentroid(const GlobalPosition& centroidComputed, const GlobalPosition& centroidExact)
{
    std::cout << std::setprecision(7) << std::scientific
              << "Centroid: " << centroidComputed << " (computed), "
                              << centroidExact << " (exact)" << "\n";
    if (Dune::FloatCmp::ne(centroidComputed, centroidExact, 1e-7))
        DUNE_THROW(Dune::Exception, "Centroid is wrong!");
}

template<class TestData>
void checkCylinderIntegration(const TestData& data, const double characteristicLength)
{
    using Point = typename TestData::Point;

    std::cout << std::endl;
    EmbeddedCoupling::CylinderIntegration<double> cylinderIntegration(characteristicLength);
    cylinderIntegration.setGeometry(data.bottomCenter, data.topCenter, data.radius, /*verbosity=*/1);
    double volumeIntegral = 0.0;
    Point centroid({0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < cylinderIntegration.size(); ++i)
    {
        volumeIntegral += cylinderIntegration.integrationElement(i);
        centroid.axpy(cylinderIntegration.integrationElement(i), cylinderIntegration.integrationPoint(i));
    }
    centroid /= volumeIntegral;

    checkVolumeIntegral(volumeIntegral, data.exactVolumeIntegral);
    checkCentroid(centroid, data.exactCentroid);
}

template<class TestData, class EllipseData>
void checkEllipticCylinderIntegration(const TestData& data, const EllipseData& ell, const double characteristicLength)
{
    using Point = typename TestData::Point;

    std::cout << std::endl;
    EmbeddedCoupling::EllipticCylinderIntegration<double> ellCylIntegration(characteristicLength);
    ellCylIntegration.setGeometry(data.bottomCenter, data.topCenter, ell.firstAxis, ell.secondAxis, /*verbosity=*/1);
    double volumeIntegral = 0.0;
    Point centroid({0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < ellCylIntegration.size(); ++i)
    {
        volumeIntegral += ellCylIntegration.integrationElement(i);
        centroid.axpy(ellCylIntegration.integrationElement(i), ellCylIntegration.integrationPoint(i));
    }
    centroid /= volumeIntegral;

    checkVolumeIntegral(volumeIntegral, data.exactVolumeIntegral);
    checkCentroid(centroid, data.exactCentroid);

    // check top cap ellipse
    std::cout << std::endl;
    EmbeddedCoupling::EllipseIntegration<double> ellIntegration(characteristicLength);
    ellIntegration.setGeometry(data.topCenter, ell.firstAxis, ell.secondAxis, /*verbosity=*/1);
    double areaIntegral = 0.0;
    Point topCentroid({0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < ellIntegration.size(); ++i)
    {
        areaIntegral += ellIntegration.integrationElement(i);
        topCentroid.axpy(ellIntegration.integrationElement(i), ellIntegration.integrationPoint(i));
    }
    topCentroid /= areaIntegral;

    checkAreaIntegral(areaIntegral, ell.firstAxis.two_norm()*ell.secondAxis.two_norm()*M_PI);
    checkCentroid(topCentroid, data.topCenter);
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using Point = Dune::FieldVector<double, 3>;

    struct TestData {
        using Point = Dune::FieldVector<double, 3>;
        Point bottomCenter, topCenter;
        double radius;
        double exactVolumeIntegral;
        Point exactCentroid;
    };

    TestData data;
    data.bottomCenter = Point({1.2, 2.3, 4.5});
    data.topCenter = Point({8.2, 4.3, 6.5});
    data.radius = 4.0;
    data.exactVolumeIntegral = M_PI*data.radius*data.radius*(data.topCenter-data.bottomCenter).two_norm();
    data.exactCentroid = data.bottomCenter; data.exactCentroid.axpy(0.5, data.topCenter-data.bottomCenter);

    /////////////////////////////////////////////////////////////////////
    // Determine cylinder volume and centroid by integration
    ////////////////////////////////////////////////////////////////////
    checkCylinderIntegration(data, 0.05);

    /////////////////////////////////////////////////////////////////////
    // Determine elliptic cylinder volume and centroid by integration
    // Construct an perfect cylinder (a == b == radius and orthogonal caps)
    ////////////////////////////////////////////////////////////////////
    struct EllipseData { Point firstAxis, secondAxis; };
    auto topNormal = data.topCenter-data.bottomCenter; topNormal /= topNormal.two_norm();
    auto firstAxis = crossProduct(topNormal, Point({1.0, 0.0, 0.0}));
    firstAxis /= firstAxis.two_norm();
    auto secondAxis = crossProduct(topNormal, firstAxis);
    // scale both axis with the same radius -> circle
    firstAxis *= data.radius; secondAxis *= data.radius;
    checkEllipticCylinderIntegration(data, EllipseData{firstAxis, secondAxis}, 0.05);

    /////////////////////////////////////////////////////////////////////
    // Determine elliptic cylinder volume and centroid by integration
    // Construct an elliptic cylinder (a == 4*b and orthogonal caps)
    ////////////////////////////////////////////////////////////////////

    firstAxis *= 2; secondAxis *= 0.5;
    checkEllipticCylinderIntegration(data, EllipseData{firstAxis, secondAxis}, 0.025);

    /////////////////////////////////////////////////////////////////////
    // Determine elliptic cylinder volume and centroid by integration
    // Construct an elliptic cylinder (a == 2*b and tilted caps)
    ////////////////////////////////////////////////////////////////////

    // rotate about the first major axis
    const auto rotationAngle = M_PI*0.25;
    auto rotationAxis = firstAxis; rotationAxis /= rotationAxis.two_norm();
    const auto rotate = make3DTransformation(1.0, Point(0.0), rotationAxis, rotationAngle, /*verbose*/false);
    secondAxis = rotate(secondAxis);
    // scale second axis so that the volume stays the same as for the orthonal cap case
    secondAxis /= std::cos(rotationAngle);
    checkEllipticCylinderIntegration(data, EllipseData{firstAxis, secondAxis}, 0.025);

    return 0;
}
