// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A simple Navier-Stokes test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

public:
    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        // determine whether to simulate a vertical or horizontal flow configuration
        verticalFlow_ = problemName_.find("vertical") != std::string::npos;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

   /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }
    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        const auto& globalPos = scvf.dofPosition();

        if (verticalFlow_)
        {
            // inflow
            if(onUpperBoundary_(globalPos))
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
                values.setDirichlet(Indices::conti0EqIdx + 1);
            }

            // left/right wall
            if (onRightBoundary_(globalPos) || (onLeftBoundary_(globalPos)))
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
                values.setNeumann(Indices::conti0EqIdx);
                values.setNeumann(Indices::conti0EqIdx + 1);
            }

        }
        else // horizontal flow
        {
            if (onLeftBoundary_(globalPos))
            {
                values.setDirichlet(Indices::conti0EqIdx + 1);
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
            }
            else if (onRightBoundary_(globalPos))
            {
                values.setDirichlet(Indices::pressureIdx);
                values.setOutflow(Indices::conti0EqIdx + 1);
            }
            else
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
                values.setNeumann(Indices::conti0EqIdx);
                values.setNeumann(Indices::conti0EqIdx + 1);
            }
        }

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::conti0EqIdx + 1);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(globalPos);

        if (verticalFlow_)
        {
            // Check if this a pure diffusion problem.
            static const bool isDiffusionProblem = problemName_.find("diffusion") != std::string::npos;

            Scalar topMoleFraction = 1e-3;

            if (isDiffusionProblem)
            {
                // For the diffusion problem, change the top mole fraction after some time
                // in order to revert the concentration gradient.
                if (time() >= 1e10)
                    topMoleFraction = 0.0;
            }
            else // advection problem
            {
                // reverse the flow direction after some time for the advection problem
                if (time() >= 3e5)
                    values[Indices::velocityYIdx] *= -1.0;
            }

            if(globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                values[Indices::conti0EqIdx + 1] = topMoleFraction;
        }
        else // horizontal flow
        {
            static const Scalar inletMoleFraction = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletMoleFraction");
            if(globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
                values[Indices::conti0EqIdx + 1] = inletMoleFraction;
        }


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

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto tmp = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::conti0EqIdx] = tmp[0];
            values[Indices::conti0EqIdx + 1] = tmp[1];
        }
        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

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
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1e5;

        static const Scalar vMax = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Velocity", 0.0);

        auto parabolicProfile = [&](const GlobalPosition& globalPos, int coord)
        {
            return vMax * (globalPos[coord] - this->gridGeometry().bBoxMin()[coord])
                        * (this->gridGeometry().bBoxMax()[coord] - globalPos[coord])
                        / (0.25 * (this->gridGeometry().bBoxMax()[coord] - this->gridGeometry().bBoxMin()[coord])
                        * (this->gridGeometry().bBoxMax()[coord] - this->gridGeometry().bBoxMin()[coord]));
        };

        if (verticalFlow_)
            values[Indices::velocityYIdx] = parabolicProfile(globalPos, 0);
        else // horizontal flow
            values[Indices::velocityXIdx] = parabolicProfile(globalPos, 1);

        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
     *        for the Beavers-Joseph-Saffman boundary condition.
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element, scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
     *        Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return 1.0;
    }

    /*!
     * \brief Sets the time loop pointer.
     */
    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    /*!
     * \brief Returns the time.
     */
    Scalar time() const
    { return timeLoop_->time(); }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    bool verticalFlow_;
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManager_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
