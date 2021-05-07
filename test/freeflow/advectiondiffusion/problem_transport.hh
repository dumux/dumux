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
/**
 * \file
 * \ingroup AdvectionDiffusionTests
 * \brief Definition of a problem for the advection diffusion transport problem:
 */

#ifndef DUMUX_ADVECTION_DIFFUSION_TEST_PROBLEM_HH
#define DUMUX_ADVECTION_DIFFUSION_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup AdvectionDiffusionTests
 *
 * \brief Transport in freeflow via the advection diffusion equation
 */
template <class TypeTag>
class AdvectionDiffusionTest : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int numFluidComps = FluidSystem::numComponents;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

public:
    AdvectionDiffusionTest(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // stating in the console whether mole or mass fractions are used
        if constexpr (useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
//     Scalar temperature() const
//     { return 273.15 + 10; } // 10C

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if(globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }


    /*!
     * \brief Evaluates the boundary conditions for a neumann
     *        boundary segment (implementation for the tpfa scheme).
     */
     PrimaryVariables neumann(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const ElementFluxVariablesCache& elemFluxVarsCache,
                              const SubControlVolumeFace& scvf) const
     {
        PrimaryVariables values(0.0);

        // no-flow everywhere except at the right boundary
        const auto xMax = this->gridGeometry().bBoxMax()[0];
        const auto& ipGlobal = scvf.ipGlobal();
        if (ipGlobal[0] < xMax - eps_)
            return values;

        // outflow condition
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const auto volFlux = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf);
        for (int compIdx = 0; compIdx < numFluidComps; ++compIdx)
        {
            values[compIdx] += volFlux / scvf.area() * ( useMoles ? (volVars.moleFraction(0, compIdx))  * volVars.molarDensity()
                                                                  : (volVars.massFraction(0, compIdx)) * volVars.density()); }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);

        if ( globalPos[0] < 0.0 + eps_
            && (globalPos[1] < 0.6 + eps_ &&
                globalPos[1] > 0.4 - eps_ ) )
        {
            if constexpr (useMoles)
                initialValues = 1e-5;
            else
                initialValues = 1e-5 * FluidSystem::molarMass(0)
                              / this->spatialParams().fluidMolarMass(globalPos);
        }
        return initialValues;
    }

    // \}

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_ = 1e-6;

    // material properties of the porous medium
    std::shared_ptr<SpatialParams> spatialParams_;
};

} // end namespace Dumux

#endif
