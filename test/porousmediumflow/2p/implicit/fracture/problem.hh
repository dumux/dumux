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
/*!
 * \file
 * \ingroup TwoPTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *        The fracture is a 2D network embedded in 3D.
 */
#ifndef DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH
#define DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include "spatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class FractureProblem;

namespace Properties
{
NEW_TYPE_TAG(FractureTypeTag, INHERITS_FROM(TwoP, FractureSpatialParams));
NEW_TYPE_TAG(FractureBoxTypeTag, INHERITS_FROM(BoxModel, FractureTypeTag));
NEW_TYPE_TAG(FractureCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, FractureTypeTag));
NEW_TYPE_TAG(FractureCCMpfaTypeTag, INHERITS_FROM(CCMpfaModel, FractureTypeTag));

// set the grid property
SET_TYPE_PROP(FractureTypeTag, Grid, Dune::FoamGrid<2, 3>);

// Set the problem property
SET_TYPE_PROP(FractureTypeTag, Problem, Dumux::FractureProblem<TypeTag>);

// Set the fluid system
SET_PROP(FractureTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Use global caching
SET_BOOL_PROP(FractureTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(FractureTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(FractureTypeTag, EnableGridFluxVariablesCache, true);

// permeablility is solution-independent
SET_BOOL_PROP(FractureTypeTag, SolutionDependentAdvection, false);
}

/*!
 * \ingroup TwoPTests
 * \brief DNAPL transport through a fracture network (2d in 3d).
 */
template <class TypeTag>
class FractureProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    enum {

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
  FractureProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
  : ParentType(fvGridGeometry) {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // // Calculate storage terms
        // PrimaryVariables storage;
        // this->model().globalStorage(storage);
        //
        // // Write mass balance information for rank 0
        // if (this->gridView().comm().rank() == 0) {
        //     std::cout<<"Storage: " << storage << std::endl;
        // }
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position where to set the BC types
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();
        if (onInlet_(globalPos))
            values.setAllNeumann();
        if (globalPos[2] > 1.0 - eps_ || globalPos[2] < eps_)
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        const auto depth = this->fvGridGeometry().bBoxMax()[dimWorld-1] - globalPos[dimWorld-1];
        const auto g = this->gravityAtPos(globalPos)[dimWorld-1];

        PrimaryVariables values;
        values[pwIdx] = 1e5 + 1000*g*depth;
        values[snIdx] = 0.0;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos)) {
            values[contiNEqIdx] = -0.04; // kg / (m * s)
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return dirichletAtPos(globalPos); }
    // \}

private:

    bool onInlet_(const GlobalPosition &globalPos) const
    { return globalPos[0] < eps_ && globalPos[1] > -0.5 - eps_; }

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} //end namespace Dumux

#endif
