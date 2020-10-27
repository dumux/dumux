// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup RichardsNCTests
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */

#ifndef DUMUX_RICHARDS_NC_WELL_TRACER_PROBLEM_HH
#define DUMUX_RICHARDS_NC_WELL_TRACER_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/richardsnc/model.hh>

#include "spatialparams.hh"

namespace Dumux {

/*!
 * \ingroup RichardsNCTests
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */
template <class TypeTag>
class RichardsWellTracerProblem;


// Specify the properties for the lens problem
namespace Properties {
// Create new type tags
namespace TTag {
struct RichardsWellTracer { using InheritsFrom = std::tuple<RichardsNC>; };
struct RichardsWellTracerBox { using InheritsFrom = std::tuple<RichardsWellTracer, BoxModel>; };
struct RichardsWellTracerCC { using InheritsFrom = std::tuple<RichardsWellTracer, CCTpfaModel>; };
} // end namespace TTag

// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsWellTracer> { using type = Dune::YaspGrid<2>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsWellTracer> { using type = RichardsWellTracerProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsWellTracer>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsWellTracerSpatialParams<GridGeometry, Scalar>;
};

// Set the physical problem to be solved
template<class TypeTag>
struct PointSource<TypeTag, TTag::RichardsWellTracer> { using type = SolDependentPointSource<TypeTag>; };
} // end namespace Properties

/*!
 * \ingroup RichardsNCTests
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards model.
 *
 * The domain is box shaped. Left and right boundaries are Dirichlet
 * boundaries with fixed water pressure (hydrostatic, gradient from right to left),
 * bottom boundary is closed (Neumann 0 boundary), the top boundary
 * (Neumann 0 boundary) is also closed. Water is extracted at a point in
 * the middle of the domain.
 * This problem is very similar to the LensProblem
 * which uses the TwoPBoxModel, with the main difference being that
 * the domain is initially fully saturated by gas instead of water and
 * water instead of a %DNAPL infiltrates from the top.
 *
 * This problem uses the \ref RichardsNCModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxrichardsnc -ParameterFile test_boxrichardsnc.input -TimeManager.TEnd 10000000</tt>
 * <tt>./test_ccrichardsnc -ParameterFile test_ccrichardsnc.input -TimeManager.TEnd 10000000</tt>
 *
 * where the initial time step is 100 seconds, and the end of the
 * simulation time is 10,000,000 seconds (115.7 days)
 */
template <class TypeTag>
class RichardsWellTracerProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    enum {
        pressureIdx = Indices::pressureIdx,
        compIdx = Indices::compMainIdx + 1,
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        dimWorld = GridView::dimensionworld
    };
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    RichardsWellTracerProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        contaminantMoleFraction_ = getParam<Scalar>("Problem.ContaminantMoleFraction");
        pumpRate_ = getParam<Scalar>("Problem.PumpRate"); // in kg/s

        // for initial conditions
        const Scalar sw = 0.4; // start with 40% saturation on top
        pcTop_ = this->spatialParams().fluidMatrixInteractionAtPos(this->gridGeometry().bBoxMax()).pc(sw);

        // for post time step mass balance
        accumulatedSource_ = 0.0;
    }

    void printTracerMass(const SolutionVector& curSol,
                         const GridVariables& gridVariables,
                         const Scalar timeStepSize)

    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;

        // bulk elements
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                tracerMass += volVars.massFraction(liquidPhaseIdx, compIdx)*volVars.density(liquidPhaseIdx)
                              * scv.volume() * volVars.saturation(liquidPhaseIdx) * volVars.porosity() * volVars.extrusionFactor();

                accumulatedSource_ += this->scvPointSources(element, fvGeometry, elemVolVars, scv)[compIdx]
                                       * scv.volume() * volVars.extrusionFactor()
                                       * FluidSystem::molarMass(compIdx)
                                       * timeStepSize;
            }
        }

        std::cout << "\033[1;33m" << "The domain contains " << tracerMass*1e9 << " µg tracer, "
                  <<  accumulatedSource_*1e9 << " µg ("<< int(std::round(-accumulatedSource_/(tracerMass - accumulatedSource_)*100))
                  <<"%) was already extracted (balanced: "
                  <<  (tracerMass - accumulatedSource_)*1e9 << " µg)\033[0m" << '\n';

    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; }; // -> 10°C

    /*!
     * \brief Returns the reference pressure [Pa] of the nonwetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonwettingReferencePressure() const
    { return 1.0e5; };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * \param globalPos The position for which the Neumann value is set
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Applies a vector of point sources which are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute rate values in units
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        auto globalPos = this->gridGeometry().bBoxMax()-this->gridGeometry().bBoxMin();
        globalPos *= 0.5;
        //! Add point source in middle of domain
        pointSources.emplace_back(globalPos,
             [this](const Problem &problem,
                    const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const ElementVolumeVariables &elemVolVars,
                    const SubControlVolume &scv)
                    {
                        const auto& volVars = elemVolVars[scv];
                        //! convert pump rate from kg/s to mol/s
                        //! We assume we can't keep up the pump rate if the saturation sinks
                        const Scalar value = pumpRate_*volVars.molarDensity(liquidPhaseIdx)/volVars.density(liquidPhaseIdx)*volVars.saturation(liquidPhaseIdx);
                        return PrimaryVariables({-value, -value*volVars.moleFraction(liquidPhaseIdx, compIdx)});
                    });
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); };

    // \}

private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        const auto xTracer = [&,this]()
        {
            const GlobalPosition contaminationPos({0.2*this->gridGeometry().bBoxMax()[0], 0.5*this->gridGeometry().bBoxMax()[1]});
            if ((globalPos - contaminationPos).two_norm() < 0.1*(this->gridGeometry().bBoxMax()-this->gridGeometry().bBoxMin()).two_norm() + eps_)
                return contaminantMoleFraction_;
            else
                return 0.0;
        }();

        PrimaryVariables values(0.0);
        //! Hydrostatic pressure profile
        values[pressureIdx] = (nonwettingReferencePressure() - pcTop_)
                               - 9.81*1000*(globalPos[dimWorld-1] - this->gridGeometry().bBoxMax()[dimWorld-1]);
        values[compIdx] = xTracer;
        return values;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
    Scalar contaminantMoleFraction_;
    Scalar pumpRate_;
    Scalar pcTop_;
    Scalar accumulatedSource_;
};

} // end namespace Dumux

#endif
