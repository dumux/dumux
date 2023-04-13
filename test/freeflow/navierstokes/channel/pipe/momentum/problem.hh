// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_FREEFLOW_PIPE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_PIPE_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {
/*!
 * \ingroup NavierStokesTests
 * \brief Freeflow problem for pipe flow
 * Simulation of a radially-symmetric pipe flow with circular cross-section
 */
template <class TypeTag, class BaseProblem>
class FreeFlowPipeProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using BoundaryTypes = typename ParentType::BoundaryTypes;

public:
    using Indices = typename ModelTraits::Indices;

    FreeFlowPipeProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        meanInletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.MeanInletVelocity");
        const Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        mu_ = nu*rho_;

        pipeRadius_ = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        pipeLength_ = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        eps_ = 1e-7*pipeRadius_;

        std::cout << "-- Reynolds number: " << 2*pipeRadius_*meanInletVelocity_/nu << std::endl;
    }

    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (onLowerBoundary_(globalPos) || onOuterBoundary_(globalPos))
                values.setAllDirichlet();
            else if (onInnerBoundary_(globalPos))
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setNeumann(Indices::momentumYBalanceIdx);
            }
            else
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();

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
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        if constexpr (ParentType::isMomentumProblem())
        {
            if (onInnerBoundary_(globalPos))
                values = 0.0; // zero shear stress at symmetry axis
            else if (onUpperBoundary_(globalPos))
                values.axpy(analyticalPressure(globalPos), scvf.unitOuterNormal());
        }
        else
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }

    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values(0.0);

        // paraboloid velocity profile
        if constexpr (ParentType::isMomentumProblem())
        {
            const auto r = globalPos[0] - this->gridGeometry().bBoxMin()[0];
            values[Indices::velocityXIdx] = 0.0;
            values[Indices::velocityYIdx] = 2.0*meanInletVelocity_*(1.0 - r*r/(pipeRadius_*pipeRadius_));
        }
        else
            values[Indices::pressureIdx] = analyticalPressure(globalPos);

        return values;
    }

    Scalar densityAtPos(const GlobalPosition& globalPos) const
    { return rho_; }

    Scalar effectiveViscosityAtPos(const GlobalPosition& globalPos) const
    { return mu_; }

    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return analyticalPressure(globalPos); }

    Scalar analyticalPressure(const GlobalPosition& globalPos) const
    {
        const auto y = globalPos[1];
        return (pipeLength_-y)*meanInletVelocity_*8.0*mu_/(pipeRadius_*pipeRadius_);
    }

    Scalar analyticalMassFlux() const
    {
        return meanInletVelocity_ * M_PI * pipeRadius_*pipeRadius_ * rho_;
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
    Scalar initialPressure_;
    Scalar meanInletVelocity_;
    Scalar mu_;
    Scalar rho_;
    Scalar pipeRadius_, pipeLength_;
    Scalar eps_;
};

} // end namespace Dumux

#endif
