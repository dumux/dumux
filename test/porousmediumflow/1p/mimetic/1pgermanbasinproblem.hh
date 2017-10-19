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
 *
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_1PGERMANBASIN_PROBLEM_HH
#define DUMUX_1PGERMANBASIN_PROBLEM_HH

#include <dune/alugrid/grid.hh>

#if PROBLEM==1
#include <dumux/implicit/staggered/properties.hh>
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#elif PROBLEM==2
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#endif
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/unit.hh>

#include "1pgermanbasinspatialparams.hh"

#include <dune/geometry/quadraturerules.hh>

namespace Dumux
{
template <class TypeTag>
class OnePGermanBasinProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<OnePGermanBasinProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
#if PROBLEM==1
NEW_TYPE_TAG(OnePGermanBasinProblem, INHERITS_FROM(StaggeredModel, OnePMimetic));
#elif PROBLEM==2
NEW_TYPE_TAG(OnePGermanBasinProblem, INHERITS_FROM(CCMpfaModel, OneP));
#endif

SET_PROP(OnePGermanBasinProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the grid type
//SET_TYPE_PROP(OnePGermanBasinProblem, Grid, Dune::YaspGrid<2>);
//SET_TYPE_PROP(OnePGermanBasinProblem, Grid, Dune::UGGrid<2>);
SET_TYPE_PROP(OnePGermanBasinProblem, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePGermanBasinProblem, Problem, Dumux::OnePGermanBasinProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePGermanBasinProblem, SpatialParams, Dumux::OnePGermanBasinSpatialParams<TypeTag> );

// Linear solver settings
//SET_TYPE_PROP(OnePGermanBasinProblem, LinearSolver, Dumux::UMFPackBackend<TypeTag> );
//SET_TYPE_PROP(OnePGermanBasinProblem, LinearSolver, Dumux::AMGBackend<TypeTag> );
SET_TYPE_PROP(OnePGermanBasinProblem, LinearSolver, Dumux::ILUnBiCGSTABBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(OnePGermanBasinProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(OnePGermanBasinProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(OnePGermanBasinProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(OnePGermanBasinProblem, EnableGlobalVolumeVariablesCache, true);
}

template<class TypeTag>
class OnePGermanBasinProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;


    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePGermanBasinProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0);
    }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;

        const GlobalPosition& normal = scvf.unitOuterNormal();

        if(std::abs(normal[dimWorld-1]) > 1.0e-2)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param intersection The intersection for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const GlobalPosition& normal = scvf.unitOuterNormal();
        if(normal[dimWorld-1] > 1.0e-5)
        {
            PrimaryVariables values(281.15);
            return values;
        }
        else
        {
            PrimaryVariables values(423.15);
            return values;
        }
    }


    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        Scalar z = globalPos[dim-1];
        PrimaryVariables priVars(281.15*(z+6071.39)/(11496.0+6071.39) - 423.15*(z-11496.0)/(11496.0+6071.39));

        return priVars;
    }

    bool shouldWriteOutput() const
    {
        return
            this->timeManager().willBeFinished();
    }

    // \}

private:
    std::string name_;
};
} //end namespace

#endif
