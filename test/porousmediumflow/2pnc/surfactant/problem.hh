#ifndef DUMUX_TEST_2P3C_SURFACTANT_PROBLEM_HH
#define DUMUX_TEST_2P3C_SURFACTANT_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag>
class TestSurfactantProblem
: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    static constexpr int dim = GridView::dimension;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

public:
    TestSurfactantProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        const auto maxSurfC = getParam<Scalar>("Problem.InjectionSurfactantConcentration");
        injectionFluidState_.setMoleFraction(0, FluidSystem::surfactantCompIdx, maxSurfC);
        this->spatialParams().setMaxSurfactantConcentration(maxSurfC);

        initialPressure_ = getParam<Scalar>("Problem.InitialPressure");
        initialSw_ = getParam<Scalar>("Problem.InitialSw");
        productionWellPressure_ = getParam<Scalar>("Problem.ProductionWellPressure");
        injectionWellPressure_ = getParam<Scalar>("Problem.InjectionWellPressure");

        name_ = getParam<std::string>("Problem.Name");
    }

    const std::string& name() const
    { return name_; }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const Scalar pw = volVars.pressure(0);
        const auto K = this->spatialParams().permeability(element);

        const auto middleX = 0.5*(this->gridGeometry().bBoxMin()[0] + this->gridGeometry().bBoxMax()[0]);
        if (scvf.ipGlobal()[0] < middleX)
        {
            using std::max;

            const Scalar dp0dn = max(0.0, (injectionWellPressure_ - volVars.pressure(0)));

            const Scalar injectionViscosity = FluidSystem::viscosity(injectionFluidState_, 0);

            // m3 / (m2 s) or kg / (m2 s)
            const auto flux = K * (1 / injectionViscosity) * dp0dn; // volume per second
            const auto moleflux = flux * FluidSystem::molarDensity(injectionFluidState_, 0);

            for (int j = 0; j < FluidSystem::numComponents; j++)
                values[j] = -moleflux * injectionFluidState_.moleFraction(0, j);

        }

        else
        {
            using std::min;
            const auto dp0dn = min(0.0, productionWellPressure_ - pw);
            const auto dp1dn = dp0dn;

            const Scalar fluxOilphase = K * volVars.mobility(1) * dp1dn; // volume per second
            const Scalar fluxWaterphase = K * volVars.mobility(0) * dp0dn; // volume per second

            values[FluidSystem::oilCompIdx] = -volVars.molarDensity(1) * fluxOilphase;
            values[FluidSystem::waterCompIdx] = -volVars.molarDensity(0) * fluxWaterphase;

            const Scalar f = volVars.moleFraction(0, FluidSystem::surfactantCompIdx);
            values[FluidSystem::surfactantCompIdx] = - volVars.molarDensity(0) * volVars.mobility(0) * K * dp0dn * f;
            values[FluidSystem::waterCompIdx] -= values[FluidSystem::surfactantCompIdx];
        }

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1, "must have p0s1 formulation");

        PrimaryVariables values(0);

        values[0] = initialPressure_;
        values[1] = 1.0 - initialSw_;
        values.setState(ModelTraits::Indices::bothPhases);
        values[FluidSystem::surfactantCompIdx] = 0;

        return values;
    }

private:
    Scalar productionWellPressure_;
    Scalar injectionWellPressure_;
    Scalar initialPressure_;
    Scalar initialSw_;

    FluidState injectionFluidState_;
    std::string name_;
};

} // end namespace Dumux

#endif
