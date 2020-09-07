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

#ifndef DUMUX_TEST_FREEFLOW_PIPE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_PIPE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

namespace Dumux {

/*!
 * \brief Freeflow problem for pipe flow
 * Simulation of a radially-symmetric pipe flow with circular cross-section
 */
template <class TypeTag>
class FreeFlowPipeProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = NavierStokesBoundaryTypes<PrimaryVariables::size()>;

public:
    FreeFlowPipeProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)

    {
        name_ = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        meanInletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.MeanInletVelocity");
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity")*getParam<Scalar>("Component.LiquidDensity");

        pipeRadius_ = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        pipeLength_ = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        eps_ = 1e-7*pipeRadius_;

        std::cout << "-- Reynolds number: " << 2*pipeRadius_*meanInletVelocity_/getParam<Scalar>("Component.LiquidKinematicViscosity") << std::endl;
    }

    const std::string& name() const
    {
        return name_;
    }

    Scalar temperature() const
    { return 293.15; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        const auto& globalPos = scvf.dofPosition();

        // inlet
        if (onLowerBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        //  outlet
        else if (onUpperBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
        }
        // pipe centerline
        else if (onInnerBoundary_(globalPos))
        {
            values.setAllSymmetry();
        }
        // pipe wall
        else if (onOuterBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        // paraboloid velocity profile
        const auto r = globalPos[0] - this->gridGeometry().bBoxMin()[0];
        const auto y = globalPos[1] - this->gridGeometry().bBoxMin()[1];
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 2.0*meanInletVelocity_*(1.0 - r*r/(pipeRadius_*pipeRadius_));
        values[Indices::pressureIdx] = (pipeLength_-y)*meanInletVelocity_*8.0*mu_/(pipeRadius_*pipeRadius_);

        return values;
    }

private:
    bool onInnerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onOuterBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    std::string name_;
    Scalar meanInletVelocity_;
    Scalar mu_;
    Scalar pipeRadius_, pipeLength_;
    Scalar eps_;
};

} // end namespace Dumux

#endif
