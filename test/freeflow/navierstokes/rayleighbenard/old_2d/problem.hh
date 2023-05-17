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
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * Flow from left to right in a two-dimensional channel is considered. At the inlet (left),
 * fixed values for velocity are set, while at the outlet (right), a fixed pressure
 * boundary condition is used. The channel is confined by solid walls at the top and bottom
 * of the domain which corresponds to no-slip/no-flow conditions.
 * For the non-isothermal test, water of increased temperature is injected at the inlet
 * while the walls are fully isolating.
 */
template <class TypeTag>
class RayleighBenardTestProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    RayleighBenardTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        temperatureBot_ = getParam<Scalar>("Problem.TemperatureBottom");
        temperatureTop_ = getParam<Scalar>("Problem.TemperatureTop");
        thermalExpansion_ = getParam<Scalar>("Problem.ThermalExpansion");
    }

    using ParentType::source;
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementFaceVariables& elemFaceVars,
                       const SubControlVolumeFace& scvf) const
    {
        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
#if NONISOTHERMAL
        source[Indices::velocityYIdx] = -this->gravity()[1] * thermalExpansion_ * volVars.priVar(Indices::temperatureIdx);
#endif
        return source;
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);
#if NONISOTHERMAL
        values.setNeumann(Indices::energyEqIdx);
        if (isPlate_(globalPos))
            values.setDirichlet(Indices::energyEqIdx);
#endif

        return values;
    }

    /*!
      * \brief Evaluates the boundary conditions for a Dirichlet control volume.
      *
      * \param globalPos The center of the finite volume which ought to be set.
      */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;
#if NONISOTHERMAL
        if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_)
            values[Indices::temperatureIdx] = temperatureBot_;
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values[Indices::temperatureIdx] = temperatureTop_;
#endif
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
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        return values;
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     * \param time A parameter for consistent signatures. It is ignored here as this analytical solution is for a stationary version of the test only.
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        PrimaryVariables values(0.0);

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

        values[Indices::pressureIdx] = 1.0e5;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

#if NONISOTHERMAL
        const static auto diag = this->gridGeometry().bBoxMax() - this->gridGeometry().bBoxMin();
        values[Indices::temperatureIdx] = temperatureBot_ + (temperatureTop_ - temperatureBot_)*globalPos[1]/diag[1];
        if ((diag/2.0 + this->gridGeometry().bBoxMin() - globalPos).two_norm() < std::abs(diag[1])*0.1)
            values[Indices::temperatureIdx] = (temperatureBot_ + temperatureTop_)/2.0;
#endif

        return values;
    }

    // \}
    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

    bool hasAnalyticalSolution()
    {
        return false;
    }

private:
    bool isPlate_(const GlobalPosition& globalPos) const
    {
        return (globalPos[1] < eps_ || globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_);
    }

    static constexpr Scalar eps_=1e-6;
    Scalar temperatureBot_;
    Scalar temperatureTop_;
    Scalar thermalExpansion_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
