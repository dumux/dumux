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
 * \ingroup RANSTests
 * \brief Channel flow test for the staggered grid RANS model.
 */

#ifndef DUMUX_RANS_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_RANS_CHANNEL_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/turbulence/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
// #include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>

namespace Dumux {

/*!
 * \ingroup RANSTests
 * \brief Channel flow test for the staggered grid RANS model.
 */
template <class TypeTag>
class RansTestProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using NumEqVector = typename ParentType::NumEqVector;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using BoundaryTypes = typename ParentType::BoundaryTypes;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    RansTestProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

            if (isOutlet_(globalPos))
                values.setAllNeumann();
        }
        else
        {
            values.setAllNeumann();
            if (isInlet_(globalPos))
            {
                values.setDirichlet(Indices::pressureIdx);
                values.setDirichlet(Indices::turbulentKineticEnergyIdx);
                values.setDirichlet(Indices::dissipationIdx);
            }

            if (isLowerWall_(globalPos) || isUpperWall_(globalPos)) //  All walls
                values.setWall();
        }

        return values;
    }


    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        PrimaryVariables values = initialAtPos(globalPos);

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = 0.0;
            if (isInlet_(globalPos))
                values[Indices::velocityXIdx] = inletVelocity_;
        }
        else
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;
            values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, true);
        }
        else
        {
            using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
            FluxVariables fluxVars;
            fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            if (isOutlet_(scvf.ipGlobal()))
                values = fluxVars.flux();
        }

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
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityYIdx] = 0.0;
            values[Indices::velocityXIdx] = inletVelocity_;
        }
        else
            values[Indices::pressureIdx] = outletPressure_;

        return values;
    }

    // \}

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return outletPressure_; }

    void setTime(Scalar time)
    { time_ = time; }

    Scalar time() const
    { return time_; }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isLowerWall_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool isUpperWall_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    static constexpr Scalar eps_=1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;
    Scalar time_;
};
} // end namespace Dumux

#endif
