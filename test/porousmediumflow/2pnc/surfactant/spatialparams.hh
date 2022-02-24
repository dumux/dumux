#ifndef DUMUX_TEST_2P3C_SURFACTANT_SPATIALPARAMS_HH
#define DUMUX_TEST_2P3C_SURFACTANT_SPATIALPARAMS_HH

#include <algorithm>

#include <dune/common/fmatrix.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/discretization/evalgradients.hh>

#include "materiallaw.hh"

namespace Dumux {

template<class GridGeometry, class Scalar>
class TestSurfactantSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                         TestSurfactantSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVPorousMediumFlowSpatialParamsMP<
        GridGeometry, Scalar, TestSurfactantSpatialParams<GridGeometry, Scalar>
    >;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int surfactantCompIdx = 2;

    using MaterialLaw = FluidMatrix::SurfactantPcKrSw<Scalar>;

public:
    using PermeabilityType = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    TestSurfactantSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        sMinKr_ = getParam<Scalar>("SpatialParams.Smin");
        sMaxKr_ = getParam<Scalar>("SpatialParams.Smax");
        srwSurf_ = getParam<Scalar>("SpatialParams.SrwSurf");
        sroSurf_ = getParam<Scalar>("SpatialParams.SroSurf");
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 0.18;
    }

    PermeabilityType permeability(const Element& element) const
    {
        return 1e-13;
    }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 300.0;
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability(element);
    }

    void setMaxSurfactantConcentration (Scalar c)
    {
        maxConcentration_ = c;
    }

    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        const auto& priVars = elemSol[scv.localDofIndex()];
        const Scalar f = priVars[surfactantCompIdx];
        const auto& gg = this->gridGeometry();

        const auto gradP = evalGradients(
            element, element.geometry(), gg, elemSol, scv.dofPosition(), true
        )[0];

        const Scalar gradPNorm = gradP.two_norm();

        using std::cos;
        const Scalar sigma =
            0.5 * 0.03 * (1 - 1e-3) * cos(f * M_PI / maxConcentration_) +
            0.5 * 0.03 * (1 + 1e-3);

        const Scalar Ncv = permeability(element, scv, elemSol) * gradPNorm / sigma;
        return MaterialLaw(Ncv, sMinKr_, sMaxKr_, srwSurf_, sroSurf_);
    }


    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return 0;
    }

private:
    Scalar maxConcentration_;
    Scalar sMinKr_, sMaxKr_, srwSurf_, sroSurf_;
};

} // end namespace Dumux

#endif
