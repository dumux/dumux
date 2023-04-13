// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Pore-network sub-problem for the coupled 1p_1p free-flow/pore-network-model test
 */

#ifndef DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_PROBLEM_PNM_HH
#define DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_PROBLEM_PNM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief  Pore-network sub-problem for the coupled 1p_1p free-flow/pore-network-model test
 *         A two-dimensional pore-network region coupled to a free-flow model.
 */
template <class TypeTag>
class PNMOnePProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    template<class SpatialParams>
    PNMOnePProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<SpatialParams> spatialParams,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, spatialParams, "PNM"), couplingManager_(couplingManager)
    {
        singleThroatTest_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.SingleThroatTest", true);
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

    const std::string& name() const
    {
        static const std::string problemName = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        return problemName;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be used for
     *        which equation for a finite volume on the boundary.
     */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        if (couplingManager().isCoupled(CouplingManager::poreNetworkIndex, CouplingManager::freeFlowMassIndex, scv))
            bcTypes.setAllCouplingNeumann();

        if (singleThroatTest_ && scv.dofIndex() == 0)
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     */
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        static const Scalar pressureBottom = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureBottom", 10.0);
        return PrimaryVariables(pressureBottom);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);

        if (couplingManager().isCoupled(CouplingManager::poreNetworkIndex, CouplingManager::freeFlowMassIndex, scv))
        {
            values = couplingManager().massCouplingCondition(
                CouplingManager::poreNetworkIndex,
                CouplingManager::freeFlowMassIndex, fvGeometry,
                scv,
                elemVolVars) / this->gridGeometry().poreVolume(scv.dofIndex());
        }

        return values;
    }

    // \}

    // \}

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    bool singleThroatTest_;
};

} // end namespace Dumux

#endif
