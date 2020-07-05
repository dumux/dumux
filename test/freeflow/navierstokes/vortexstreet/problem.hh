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

#ifndef DUMUX_TEST_KARMAN_VORTEX_STREET_PROBLEM_HH
#define DUMUX_TEST_KARMAN_VORTEX_STREET_PROBLEM_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/material/components/air.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief The Karman vortex street problem for flow around a cylinder
 */
template <class TypeTag>
class KarmanVortexStreetProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = NavierStokesBoundaryTypes<PrimaryVariables::size()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    KarmanVortexStreetProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        cylinderRadius_ = getParam<Scalar>("Grid.CylinderRadius");
        cylinderPos_ = getParam<GlobalPosition>("Grid.CylinderPosition");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.0e5);

        const auto Re = getParam<Scalar>("Problem.ReynoldsNumber");
        using Air = Components::Air<Scalar>;
        const auto nu = Air::gasViscosity(temperature(), outletPressure_)/Air::gasDensity(temperature(), outletPressure_);

        inletVelocity_ = Re*nu/(2.0*cylinderRadius_);

        std::cout << "-- Inlet velocity: " << inletVelocity_ << " m/s" << "\n";
        std::cout << "-- Reynolds number: " << Re << "\n";

        const auto St = 0.198*(1.0 - 19.7/Re);
        std::cout << "-- Approximate Strouhal number: " << St << "\n";
        std::cout << "-- Frequency estimate: " << St*inletVelocity_/(2.0*cylinderRadius_) << " Hz" << "\n";
    }

    Scalar temperature() const { return 293.15; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (isInlet_(globalPos) || onCylinderBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else if (isOutlet_(globalPos))
            values.setDirichlet(Indices::pressureIdx);

        else
            values.setAllSymmetry();

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        if (isInlet_(globalPos))
        {
            PrimaryVariables values(0.0);
            values[Indices::velocityYIdx] = 0.0;
            values[Indices::velocityXIdx] = inletVelocity_;
            return values;
        }

        if (onCylinderBoundary_(globalPos))
        {
            PrimaryVariables values(0.0);
            values[Indices::velocityYIdx] = 0.0;
            values[Indices::velocityXIdx] = 0.0;
            return values;
        }

        if (isOutlet_(globalPos))
        {
            PrimaryVariables values(0.0);
            values[Indices::pressureIdx] = outletPressure_;
            return values;
        }

        DUNE_THROW(Dune::InvalidStateException, "Called problem.dirichlet() on non-dirichlet boundary!");
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = outletPressure_;
        values[Indices::velocityYIdx] = 0.0;
        values[Indices::velocityXIdx] = inletVelocity_;
        return values;
    }

private:
    bool onCylinderBoundary_(const GlobalPosition& globalPos) const
    { return (globalPos-cylinderPos_).two_norm() < cylinderRadius_*1.5; }

    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_ = 1e-12;
    Scalar inletVelocity_;
    Scalar outletPressure_;
    Scalar cylinderRadius_;
    GlobalPosition cylinderPos_;
};

} // end namespace Dumux

#endif
