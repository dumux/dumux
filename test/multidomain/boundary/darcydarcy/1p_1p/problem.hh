// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The properties for the incompressible test.
 */

#ifndef DUMUX_ONEP_SUB_TEST_PROBLEM_HH
#define DUMUX_ONEP_SUB_TEST_PROBLEM_HH

#include <dune/common/indices.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief Multidomain test problem for the incompressible one-phase model.
 *
 * Two possibilities to divide the model domain are given:
 * half: a horizontal interface splits the domain in two equally sized subdomains
 * lens: one subdomain is defined as a central lens, surrounded by the other subdomain
 */
template<class TypeTag, std::size_t tag>
class OnePTestProblem
: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr auto domainIdx = Dune::index_constant<tag>{};

public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        // set a default name for the problem
        problemName_ = getParam<std::string>("Vtk.OutputName")+ "_" + getParamFromGroup<std::string>(paramGroup, "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        const auto& globalPos = scvf.ipGlobal();

        if (globalPos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps_
            || globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        if (couplingManager_->isCoupled(domainIdx, scvf))
            values.setAllCouplingNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto bcTypes = boundaryTypes(element, scvf);
        if (bcTypes.hasCouplingNeumann())
            values[Indices::conti0EqIdx] = couplingManager_->advectiveFluxCoupling(domainIdx, element, fvGeometry, elemVolVars, scvf);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[0] = 1.0e+5*(2.0 - globalPos[dimWorld-1]);
        return values;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return dirichletAtPos(globalPos); }


private:
    std::shared_ptr<CouplingManager> couplingManager_;
    static constexpr Scalar eps_ = 1e-7;
    std::string problemName_;

};

} // end namespace Dumux

#endif
