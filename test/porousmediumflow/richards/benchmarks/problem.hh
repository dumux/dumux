// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Richards benchmarks base problem
 *
 * Infiltration benchmark:
 * Root-soil benchmark paper Schnepf et al. (case M2.1, Eq. 4) https://doi.org/10.3389/fpls.2020.00316
 * based on Vanderborght 2005 (see Fig. 4abc and Eq. 56-60) https://doi.org/10.2113/4.1.206
 *
 * Evaporation benchmark:
 * Root-soil benchmark paper Schnepf et al. (case M2.2) https://doi.org/10.3389/fpls.2020.00316
 * based on Vanderborght 2005 (see Fig. 5abcd and Eq. 39-47) https://doi.org/10.2113/4.1.206
 */

#ifndef DUMUX_RICHARDS_BECHMARK_PROBLEM_HH
#define DUMUX_RICHARDS_BECHMARK_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/io/format.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

enum class BenchmarkScenario {
    evaporation, infiltration
};

/*!
 * \ingroup RichardsTests
 * \brief Richards benchmarks base problem
 */
template <class TypeTag>
class RichardsBenchmarkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    RichardsBenchmarkProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        const auto density = Components::SimpleH2O<double>::liquidDensity(0,0);
        const auto initialHead = getParam<Scalar>("Problem.InitialHeadInCm")*0.01;
        const auto criticalHead = getParam<Scalar>("Problem.CriticalSurfaceHeadInCm")*0.01;
        initialPressure_ = 1.0e5 + initialHead*9.81*density;
        criticalSurfacePressure_ = 1.0e5 + criticalHead*9.81*density;
        const auto& fm = this->spatialParams().fluidMatrixInteractionAtPos(0.0);
        const auto criticalSaturation = fm.sw(1.0e5 - criticalSurfacePressure_);
        criticalSurfaceKrw_ = fm.krw(criticalSaturation);
        enableGravity_ = getParam<bool>("Problem.EnableGravity", true);

        const auto potentialRate = getParam<Scalar>("Problem.SurfaceFluxMilliMeterPerDay");
        potentialRate_ = density*potentialRate/(1000*86400.0); // mm/day -> kg/s/m^2
        useKrwAverage_ = getParam<bool>("Problem.UseKrwAverage", false);
        bottomDirichlet_ = getParam<bool>("Problem.BottomDirichlet", false);
        scenario_ = (potentialRate > 0) ? BenchmarkScenario::evaporation : BenchmarkScenario::infiltration;
        surfaceArea_ = (scenario_ == BenchmarkScenario::evaporation) ? 0.1*0.1 : 0.05*0.05;
    }

    // output name
    const std::string& name() const
    { return name_; }

    // reference pressure
    Scalar nonwettingReferencePressure() const
    { return 1.0e5; };

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onLowerBoundary(globalPos) && !bottomDirichlet_)
            bcTypes.setAllNeumann();
        else if (onLowerBoundary(globalPos) && bottomDirichlet_)
            bcTypes.setAllDirichlet();
        else if (onUpperBoundary(globalPos))
            bcTypes.setAllNeumann();
        else
            DUNE_THROW(Dune::InvalidStateException, "Wrong boundary?");
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    /*!
     * \brief Evaluates the initial values for a control volume
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = initialPressure_;
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     * Negative values mean influx.
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        if (onUpperBoundary(globalPos))
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];

            const auto dist = (fvGeometry.scv(scvf.insideScvIdx()).center() - globalPos).two_norm();
            const auto cellPressure = volVars.pressure(0);
            const auto density = volVars.density(0);
            const auto viscosity = volVars.viscosity(0);
            const auto relPerm = volVars.relativePermeability(0);
            const auto K = volVars.permeability();
            const auto gravity = enableGravity_ ? 9.81 : 0.0;
            const auto avgRelPerm = 0.5*(relPerm + criticalSurfaceKrw_);

            // kg/m^3 * m^2 * Pa / m / Pa / s = kg/s/m^2
            auto criticalRate = density*K/viscosity*((cellPressure - criticalSurfacePressure_)/dist - density*gravity);

            if (!std::signbit(criticalRate))
                criticalRate *= useKrwAverage_ ? avgRelPerm : relPerm;

            if (scenario_ == BenchmarkScenario::evaporation)
                values[Indices::conti0EqIdx] = std::min(potentialRate_, criticalRate);
            else
                values[Indices::conti0EqIdx] = std::max(potentialRate_, criticalRate);
        }

        // free drainage (gravitational flux) for infiltration scenario
        else if (onLowerBoundary(globalPos) && (scenario_ == BenchmarkScenario::infiltration))
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            const auto gravity = enableGravity_ ? 9.81 : 0.0;
            const auto density = volVars.density(0);
            const auto viscosity = volVars.viscosity(0);
            const auto relPerm = volVars.relativePermeability(0);
            const auto K = volVars.permeability();

            values[Indices::conti0EqIdx] = density*K*relPerm/viscosity*(density*gravity);
        }

        return values;
    }

    bool onLowerBoundary(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps_; }

    bool onUpperBoundary(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_; }

    //! compute the actual evaporation/infiltration rate
    template<class SolutionVector, class GridVariables>
    Scalar computeActualRate(const SolutionVector& sol, const GridVariables& gridVars, bool verbose = true) const
    {
        Scalar rate = 0.0;

        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVars.curGridVolVars());

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, sol);
            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary())
                    rate += this->neumann(element, fvGeometry, elemVolVars, 0.0, scvf)[0];
        }

        if (verbose)
            std::cout << Fmt::format("Actual rate: {:.5g} (mm/day)\n", rate*86400*1000/1000);

        return rate*86400*1000/1000;
    }

    /*!
     * \brief Adds Robin flux derivatives for wetting phase
     *
     * \param derivativeMatrices The matrices containing the derivatives
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param curElemVolVars The current element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The sub control volume face
     */
    template<class PartialDerivativeMatrices, class ElementVolumeVariables, class ElementFluxVariablesCache>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto& globalPos = scvf.ipGlobal();
        const auto insideFluidMatrixInteraction = this->spatialParams().fluidMatrixInteractionAtPos(globalPos);

        //material law derivatives
        const auto insideSw = insideVolVars.saturation(0);
        const auto insidePc = insideVolVars.capillaryPressure();
        const auto dsw_dpw_inside = -insideFluidMatrixInteraction.dsw_dpc(insidePc);
        const auto dkrw_dsw_inside = insideFluidMatrixInteraction.dkrw_dsw(insideSw);

        if (onUpperBoundary(globalPos))
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];

            const auto dist = (fvGeometry.scv(scvf.insideScvIdx()).center() - globalPos).two_norm();
            const auto cellPressure = volVars.pressure(0);
            const auto density = volVars.density(0);
            const auto viscosity = volVars.viscosity(0);
            const auto relPerm = useKrwAverage_ ?  volVars.relativePermeability(0)*0.5 : volVars.relativePermeability(0);
            const auto K = volVars.permeability();
            const auto gravity = enableGravity_ ? 9.81 : 0.0;

            // kg/m^3 * m^2 * Pa / m / Pa / s = kg/s/m^2
            auto criticalRate = density*K/viscosity*((cellPressure - criticalSurfacePressure_)/dist - density*gravity);

            if (!std::signbit(criticalRate))
                criticalRate *= relPerm;

            if (scenario_ == BenchmarkScenario::evaporation)
            {
                if (criticalRate <= potentialRate_)
                    derivativeMatrices[insideScvIdx][Indices::conti0EqIdx][0]
                        += (density/viscosity*K*relPerm/dist + density*K/viscosity*((cellPressure - criticalSurfacePressure_)/dist - density*gravity) *dkrw_dsw_inside*dsw_dpw_inside)*surfaceArea_;
            }
            else //in case of infiltration no relative permeability is added in this term
            {
                if (criticalRate >= potentialRate_)
                    derivativeMatrices[insideScvIdx][Indices::conti0EqIdx][0] += density*K/viscosity/dist*surfaceArea_;
            }
        }

        //free drainage (gravitational flux) for infiltration scenario
        else if (onLowerBoundary(globalPos) && (scenario_ == BenchmarkScenario::infiltration))
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            const auto gravity = enableGravity_ ? 9.81 : 0.0;
            const auto density = volVars.density(0);
            const auto relPerm = volVars.relativePermeability(0);
            const auto viscosity = volVars.viscosity(0);
            const auto K = volVars.permeability();

            derivativeMatrices[insideScvIdx][Indices::conti0EqIdx][0] += density*K*relPerm/viscosity*(density*gravity)*dkrw_dsw_inside*dsw_dpw_inside*surfaceArea_;
        }
    }

private:
    Scalar initialPressure_, criticalSurfacePressure_, potentialRate_;
    Scalar criticalSurfaceKrw_;
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
    bool enableGravity_;
    bool bottomDirichlet_;
    bool useKrwAverage_;
    BenchmarkScenario scenario_;
    Scalar surfaceArea_;
};

} // end namespace Dumux

#endif
