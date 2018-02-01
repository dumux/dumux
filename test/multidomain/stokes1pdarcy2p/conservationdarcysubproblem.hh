// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief \todo please doc me
 */
#ifndef DUMUX_CONSERVATION_DARCY_SUBPROBLEM_HH
#define DUMUX_CONSERVATION_DARCY_SUBPROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "conservationspatialparams.hh"
#include <dumux/material/fluidsystems/h2oair.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{
template <class TypeTag>
class ConservationDarcySubProblem;

namespace Properties
{
NEW_TYPE_TAG(DarcySubProblem, INHERITS_FROM(CCTpfaModel, TwoP, ConservationSpatialParams));

// Set the problem property
SET_TYPE_PROP(DarcySubProblem, Problem, ConservationDarcySubProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(DarcySubProblem, SpatialParams, ConservationSpatialParams<TypeTag>);

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(DarcySubProblem, Grid, Dune::YaspGrid<3>);
#else
SET_TYPE_PROP(DarcySubProblem, Grid, Dune::YaspGrid<2>);
#endif

SET_TYPE_PROP(DarcySubProblem, FluidSystem, FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Disable gravity
SET_BOOL_PROP(DarcySubProblem, ProblemEnableGravity, false);

// choose pn and Sw as primary variables
SET_INT_PROP(DarcySubProblem, Formulation, TwoPFormulation::pnsw);

// // the gas component balance (air) is replaced by the total mass balance
// SET_INT_PROP(DarcySubProblem, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, Indices)::contiNEqIdx);

SET_BOOL_PROP(DarcySubProblem, EnableGlobalFVGeometryCache, true);

// TODO define in global problem?
SET_BOOL_PROP(DarcySubProblem, UseMoles, false);

// Set the grid parameter group
SET_STRING_PROP(DarcySubProblem, GridParameterGroup, "DarcyGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief \todo please doc me
 */
template <class TypeTag >
class ConservationDarcySubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ConservationDarcySubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        pressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, Pressure);
        saturation_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, Saturation);
        temperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, Temperature);
        name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);

        bBoxMin_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, DarcyGrid, LowerLeft);
        bBoxMax_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, DarcyGrid, UpperRight);

        eps_ = 1e-6;

        this->boundingBoxTree(); // TODO
     }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     *
     */
    const std::string name() const
    {
        return name_+"_darcy";
    }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        ParentType::init();
    }

    // suppress output from DuMuX
    bool shouldWriteOutput() const
    {
        return false;
    }

    /*!
     * \brief Returns the source term at specific position in the domain.
     *
     * \param values The source values for the primary variables
     * \param globalPos The position
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        return values;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
    * \name Boundary conditions
    */
   // \{
   //! Specifies which kind of boundary condition should be used for
   //! which equation for a finite volume on the boundary.
   BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
   {
       BoundaryTypes bcTypes;
       bcTypes.setAllNeumann();

       double eps = 1.0e-5;
       if (scvf.center()[dim-1] < eps) // bottom boundary
       {
           bcTypes.setDirichlet(pressureIdx);
           bcTypes.setDirichlet(saturationIdx);
       }

       return bcTypes;
   }


   /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        control volume.
    *
    * \param values The dirichlet values for the primary variables
    * \param vertex The vertex (pore body) for which the condition is evaluated
    *
    * For this method, the \a values parameter stores primary variables.
    */
   PrimaryVariables dirichlet(const Element &element,
                              const SubControlVolumeFace &scvf) const
   {
       PrimaryVariables values(0.0);
       values[pressureIdx] = pressure_;
       values[saturationIdx] = saturation_;

       return values;
   }


    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The element volume variables
     * \param scvf The sub control volume face
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
   PrimaryVariables neumann(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolvars,
                            const SubControlVolumeFace& scvf) const
   {
       PrimaryVariables values(0.0);

       if (onCouplingInterface(scvf.center()))
       {
           values[pressureIdx] = couplingManager().stokesData().massCouplingCondition(scvf);
       }
       return values;
   }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = pressure_;
        values[saturationIdx] = saturation_;
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    Scalar temperature() const
    {
        return temperature_;
    }

    /*!
     * \brief Set the coupling manager
     *
     * \param couplingManager The coupling manager for the global problem
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Check if on coupling interface
     *
     * \param globalPos The global position
     *
     * Returns true if globalPos is on coupling interface
     * (here: upper boundary of Darcy domain)
     */
    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return globalPos[1] > bBoxMax_[1] - eps_; }

private:
    Scalar pressure_;
    Scalar saturation_;
    Scalar temperature_;
    std::string name_;

    Scalar eps_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif // DUMUX_CONSERVATION_DARCY_SUBPROBLEM_HH
