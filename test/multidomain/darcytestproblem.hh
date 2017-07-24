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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * @file
 * @brief  Definition of a simple Darcy problem
 */
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

// base problem
// 1p porous medium flow model

// #include <dumux/porousmediumflow/1p/implicit/properties.hh> // TODO remove?
#include <test/porousmediumflow/1p/implicit/1ptestspatialparams.hh> // TODO replace!!
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>
// #include <dumux/porenetworkflow/common/functions.hh>
// #include <dumux/io/matchinggridcreator.hh> // TODO


namespace Dumux
{
template <class TypeTag>
class DarcyTestProblem;

namespace Properties
{
NEW_TYPE_TAG(DarcyTestProblem, INHERITS_FROM(CCTpfaModel, OneP, OnePTestSpatialParams)); // TODO ??

// Set the problem property
SET_TYPE_PROP(DarcyTestProblem, Problem, Dumux::DarcyTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(DarcyTestProblem, SpatialParams, OnePTestSpatialParams<TypeTag>);

// SET_TYPE_PROP(DarcyTestProblem, GridCreator, Dumux::StructuredMicroModelGridCreator<TypeTag>);

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(DarcyTestProblem, Grid, Dune::YaspGrid<3>); //, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
// SET_TYPE_PROP(DarcyTestProblem, GridCreator, MatchingGridCreator<TypeTag, 3>);
#else
SET_TYPE_PROP(DarcyTestProblem, Grid, Dune::YaspGrid<2>); //, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
// SET_TYPE_PROP(DarcyTestProblem, GridCreator, MatchingGridCreator<TypeTag, 2>);
#endif

SET_PROP(DarcyTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(DarcyTestProblem, ProblemEnableGravity, false);

// SET_BOOL_PROP(DarcyTestProblem, NeglectPoreFlowResistance, false); // TODO

// Set the grid parameter group
SET_STRING_PROP(DarcyTestProblem, GridParameterGroup, "DarcyGrid");

// SET_TYPE_PROP(DarcyTestProblem, SinglePhaseTransmissibility, TransmissibilityBruus<TypeTag, false>); // TODO

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
NEW_PROP_TAG(StokesProblemTypeTag);
}

template <class TypeTag>
class DarcyTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        // primary variable indices
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimworld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalTypeTag;
    typedef typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager) CouplingManager;
    typedef typename GET_PROP_TYPE(GlobalTypeTag, StokesProblemTypeTag) StokesProblemTypeTag;
    typedef typename GET_PROP_TYPE(StokesProblemTypeTag, Indices) StokesIndices;

    enum { dofCodim = 0 };
//     typedef typename Dumux::Functions<TypeTag> Functions;

public:
    DarcyTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-7)
    {
        // get some parameters from the input file
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {
        return name_+"_darcy";
    }

    bool shouldWriteOutput() const //define output
    {
        // write initial conditions
        return (this->timeManager().time() < 0.0);
    }

    void preTimeStep()
    {
    }

    /*!
     * \brief Called at the end of each time step
     */
    void postTimeStep()
    {
//         std::cout << "Darcy outflux: " << Functions::boundaryFlux(*this,"max", 0) << std::endl;
    }


    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

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
        if(couplingManager().isDarcyCouplingEntity(scvf))
            bcTypes.setAllCouplingNeumann();

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
        values[pressureIdx] = 1.0e5;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables neumann(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        return PrimaryVariables(0.0);
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    PrimaryVariables source(const Element& element, const VolumeVariables& volVars)
    {
        return PrimaryVariables(0.0);
    }
    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Element &element) const
    {
        PrimaryVariables values(0.0);
            values[pressureIdx] = 1.0e5;
            return values;
    }


    // \}

   void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    //! Get the coupling manager
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar eps_;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif
