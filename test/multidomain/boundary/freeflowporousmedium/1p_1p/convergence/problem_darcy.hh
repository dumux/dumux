// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/porousmediumflow/problem.hh>

#include "testcase.hh"
#include "analyticalsolutions.hh"

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test
 */
template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class BC {
        dirichlet, neumann, mixed
    };

public:
    //! export the Indices
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    const DarcyStokesTestCase testCase,
                    const std::string& name)
    : ParentType(gridGeometry, spatialParams, "Darcy")
    , couplingManager_(couplingManager)
    , testCase_(testCase)
    {
        problemName_ = name + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        auto bc = getParamFromGroup<std::string>(this->paramGroup(), "Problem.BoundaryConditions", "Dirichlet");
        if (bc == "Dirichlet")
            boundaryConditions_ = BC::dirichlet;
        else if (bc == "Neumann")
            boundaryConditions_ = BC::neumann;
        else if (bc == "Mixed")
            boundaryConditions_ = BC::mixed;
        else
            DUNE_THROW(Dune::Exception, "Wrong BC type choose: Dirichlet, Neumann or Mixed");

        std::cout << "Porous medium domain: Using " << bc << " boundary conditions" << std::endl;
    }

    const std::string& name() const
    { return problemName_; }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;

        if (couplingManager().isCoupled(CouplingManager::porousMediumIndex, CouplingManager::freeFlowMassIndex, scvf))
            values.setAllCouplingNeumann();
        else
        {
            if (boundaryConditions_ == BC::dirichlet)
                values.setAllDirichlet();
            else if (boundaryConditions_ == BC::neumann)
                values.setAllNeumann();
            else
            {
                if (onLeftBoundary_(scvf.center()))
                    values.setAllNeumann();
                else
                    values.setAllDirichlet();
            }
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary sub-control-volume-face
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const auto p = fullAnalyticalSolution(scvf.center())[2];
        return PrimaryVariables(p);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupled(CouplingManager::porousMediumIndex, CouplingManager::freeFlowMassIndex, scvf))
            values[Indices::conti0EqIdx] = couplingManager().massCouplingCondition(
                CouplingManager::porousMediumIndex, CouplingManager::freeFlowMassIndex,
                fvGeometry, scvf, elemVolVars
            );

        else
        {
            const auto sol = fullAnalyticalSolution(scvf.center());
            const auto n = scvf.unitOuterNormal();
            auto v = n; v[0] = sol[0]; v[1] = sol[1];
            values[Indices::conti0EqIdx] = v*n;
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                return ShiueOne::darcyRHS(globalPos);
            case DarcyStokesTestCase::ShiueExampleTwo:
                return ShiueTwo::darcyRHS(globalPos);
            case DarcyStokesTestCase::Rybak:
                return Rybak::darcyRHS(globalPos);
            case DarcyStokesTestCase::Schneider:
                return Schneider::darcyRHS(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     * \param element The element
     */
    PrimaryVariables initial(const Element &element) const
    {  return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     * Returns vector with entries: (velocity-x | velocity-y | pressure)
     */
    Dune::FieldVector<Scalar, 3> fullAnalyticalSolution(const GlobalPosition& globalPos) const
    {
        using namespace Solution::DarcyStokes;
        switch (testCase_)
        {
            case DarcyStokesTestCase::ShiueExampleOne:
                return ShiueOne::darcy(globalPos);
            case DarcyStokesTestCase::ShiueExampleTwo:
                return ShiueTwo::darcy(globalPos);
            case DarcyStokesTestCase::Rybak:
                return Rybak::darcy(globalPos);
            case DarcyStokesTestCase::Schneider:
                return Schneider::darcy(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    static constexpr Scalar eps_ = 1e-7;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::string problemName_;
    DarcyStokesTestCase testCase_;
    BC boundaryConditions_;
};

} // end namespace Dumux

#endif
