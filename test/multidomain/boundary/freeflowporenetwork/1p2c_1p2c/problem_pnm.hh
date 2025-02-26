// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Pore-network sub-problem for the coupled 1p2c_1p2c free-flow/pore-network-model test
 */

#ifndef DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_ONEPTWOC_PROBLEM_PNM_HH
#define DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_ONEPTWOC_PROBLEM_PNM_HH

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
class PNMOnePNCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    template<class SpatialParams>
    PNMOnePNCProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<SpatialParams> spatialParams,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, spatialParams, "PNM"), couplingManager_(couplingManager)
    {
        initialPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialPressure", 1e5);
        initialMoleFraction_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialMoleFraction", 2e-3);
        onlyDiffusion_ = getParam<bool>("Problem.OnlyDiffusion", false);
#if !ISOTHERMAL
        initialTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialTemperature", 273.15 + 20);
#endif
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
        else if (onLowerBoundary_(scv))
        {
            if (!onlyDiffusion_ )
                bcTypes.setAllDirichlet();
            else
                bcTypes.setAllNeumann();
#if !ISOTHERMAL
            bcTypes.setDirichlet(Indices::temperatureIdx);
#endif
        }

        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     */
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        return initialAtPos(scv.dofPosition());
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
            const auto& massCouplingCondition = couplingManager().massCouplingCondition(
                CouplingManager::poreNetworkIndex, CouplingManager::freeFlowMassIndex,
                fvGeometry, scv, elemVolVars
            );
            values[Indices::conti0EqIdx] = massCouplingCondition[Indices::conti0EqIdx];
            values[Indices::conti0EqIdx + 1] = massCouplingCondition[Indices::conti0EqIdx + 1];
#if !ISOTHERMAL
            values[Indices::energyEqIdx] = couplingManager().energyCouplingCondition(
                CouplingManager::poreNetworkIndex,
                CouplingManager::freeFlowMassIndex, fvGeometry,
                scv,
                elemVolVars);
#endif
        }
        values /= this->gridGeometry().poreVolume(scv.dofIndex());

        return values;
    }

    /*!
     * \brief Evaluate the initial value for a given global position.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = initialPressure_;
        values[Indices::conti0EqIdx + 1] = initialMoleFraction_;

#if !ISOTHERMAL
        values[Indices::temperatureIdx] = initialTemperature_;
#endif
        return values;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    bool onLowerBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == 2; }

    Scalar initialPressure_;
    Scalar initialMoleFraction_;
    bool onlyDiffusion_;
#if !ISOTHERMAL
    Scalar initialTemperature_;
#endif
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
