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

#ifndef DUMUX_CONSERVATION_STOKES_SUBPROBLEM_HH
#define DUMUX_CONSERVATION_STOKES_SUBPROBLEM_HH

//#include <appl/staggeredgrid/freeflow/navierstokes/navierstokes2cni/navierstokes2cniproblem.hh> TODO - ok?!

//#include "../properties.hh" // TODO

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggerednc/model.hh>
#include <dumux/implicit/problem.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{
template <class TypeTag>
class ConservationStokesSubProblem;

namespace Properties
{
NEW_TYPE_TAG(StokesSubProblem, INHERITS_FROM(StaggeredModel, NavierStokesNCNI));

// Set the problem property
SET_TYPE_PROP(StokesSubProblem, Problem, ConservationStokesSubProblem<TypeTag>);

// TODO set in navierstokes2cnipropertydefaults, different in staggerednc/propertydefaults
SET_INT_PROP(StokesSubProblem, PhaseIdx, GET_PROP_TYPE(TypeTag, FluidSystem)::nPhaseIdx);

// Set the grid type // TODO from stokestestproblem (1p)
#if ENABLE_3D
SET_TYPE_PROP(StokesSubProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
#else
SET_TYPE_PROP(StokesSubProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
#endif

// Use complete Navier-Stokes equation
SET_BOOL_PROP(StokesSubProblem, EnableInertiaTerms, false);

// Disable gravity field
SET_BOOL_PROP(StokesSubProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(StokesSubProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(StokesSubProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(StokesSubProblem, EnableGlobalVolumeVariablesCache, true);

// Set the grid parameter group
SET_STRING_PROP(StokesSubProblem, GridParameterGroup, "StokesGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}

template <class TypeTag>
class ConservationStokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    enum { // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { // primary variable indices
            pressureIdx = Indices::pressureIdx,
            velocityXIdx = Indices::velocityXIdx,
            velocityYIdx = Indices::velocityYIdx,
            velocityZIdx = Indices::velocityZIdx,
            massOrMoleFracIdx = Indices::massOrMoleFracIdx,
            temperatureIdx = Indices::temperatureIdx
    };

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

public:
    ConservationStokesSubProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView), gridView_(gridView)
    {
        velocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, Velocity);
        pressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, Pressure);
        temperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, Temperature);
        massMoleFrac_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, MassMoleFrac);

        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[0];
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[0];
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[dim-1];
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[dim-1];
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NameFF); }

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (isOnCouplingFace(globalPos))
        {
            values.setAllDirichlet();
            values.setNeumann(velocityXIdx);
            values.setNeumann(velocityYIdx);
        }
        else
        {
            values.setOutflow(pressureIdx);
            values.setOutflow(massOrMoleFracIdx); // TODO inflow
            values.setOutflow(temperatureIdx); // TODO inflow
            if (globalPos[dim-1] < bBoxMin_[dim-1] + eps_)
            {
                values.setDirichlet(velocityXIdx); // TODO inflow
                values.setDirichlet(velocityYIdx); // TODO inflow
            }
            else
            {
                values.setDirichlet(velocityXIdx); // TODO wall
                values.setDirichlet(velocityYIdx); // TODO wall
            }
        }

        return values;
    }

    /*!
     * \brief Return Dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    BoundaryValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryValues values(0.0);
        values[pressureIdx] = pressure_;

        if (dim == 1)
            values[velocityXIdx] = velocity_;
        else // dim == 2
            values[velocityYIdx] = 4.0 * velocity_ * globalPos[0] * (bBoxMax_[0] - globalPos[0])
                   / (bBoxMax_[0] - bBoxMin_[0])
                   / (bBoxMax_[0] - bBoxMin_[0]);

        values[massOrMoleFracIdx] = massMoleFrac_;

        values[temperatureIdx] = temperature_;

        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0);

        values = dirichletAtPos(globalPos);

        return values;
    }

    //! \brief Returns whether we are on a coupling face
    const bool isOnCouplingFace(const GlobalPosition& global) const
    {
        if (dim == 1)
            return global[0] < bBoxMin_[0] + eps_;
        else // dim == 2
            return global[dim-1] < bBoxMin_[dim-1] + eps_;
    }

    /*!
     * \brief Set the coupling manager
     * \param couplingManager The coupling manager
     *
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     *
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return globalPos[1] < bBoxMin_[1] + eps_; } // onLowerBoundary

private:
    std::string name_;
    const GridView gridView_;
    static constexpr Scalar eps_ = 1e-6;

    Scalar velocity_;
    Scalar pressure_;
    Scalar temperature_;
    Scalar massMoleFrac_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif // DUMUX_CONSERVATION_STOKES_SUBPROBLEM_HH
