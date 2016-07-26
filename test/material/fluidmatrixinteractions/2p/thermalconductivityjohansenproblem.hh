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
 * \brief Simple test problem for the Johansen thermal conductivity law
 */
#ifndef DUMUX_THERMAL_CONDUCTIVITY_JOHANSEN_PROBLEM_HH
#define DUMUX_THERMAL_CONDUCTIVITY_JOHANSEN_PROBLEM_HH

#include <dumux/material/fluidsystems/h2on2.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivityjohansen.hh>

#include "thermalconductivityspatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class ThermalConductivityJohansenProblem;

namespace Properties
{
NEW_TYPE_TAG(ThermalConductivityJohansenProblem, INHERITS_FROM(BoxModel, TwoPTwoCNI, ThermalConductivitySpatialParams));

// Set the grid type
SET_TYPE_PROP(ThermalConductivityJohansenProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ThermalConductivityJohansenProblem, Problem, ThermalConductivityJohansenProblem<TypeTag>);

// Set the wetting phase
SET_TYPE_PROP(ThermalConductivityJohansenProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

// Set thermal conductivity law
SET_TYPE_PROP(ThermalConductivityJohansenProblem, ThermalConductivityModel,
              ThermalConductivityJohansen<typename GET_PROP_TYPE(TypeTag, Scalar)>);
}


/*!
 * \ingroup MaterialTestProblems
 *
 * \brief Simple test problem for the Johansen thermal conductivity law
 *
 * To run the test execute the following line in shell:
 * <tt>./test_thermalconductivityjohansen</tt>
 *
 */
template <class TypeTag >
class ThermalConductivityJohansenProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
        temperatureIdx = Indices::temperatureIdx,

        // Phase State
        wPhaseOnly = Indices::wPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ThermalConductivityJohansenProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        FluidSystem::init();
        this->spatialParams().plotMaterialLaw();
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
    const std::string name() const
    { return "test_thermalconductivityjohansen"; }


    //! \copydoc ImplicitProblem::sourceAtPos()
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{


    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllDirichlet();
    }


    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }


    //! \copydoc ImplicitProblem::neumann()
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    //! \copydoc ImplicitProblem::initialAtPos()
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }


    //! \copydoc InjectionProblem::initialPhasePresence()
    int initialPhasePresence(const Vertex &vertex,
                             int &vIdxGlobal,
                             const GlobalPosition &globalPos) const
    {
        return wPhaseOnly;
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + globalPos[1]*densityW*9.81;
        values[switchIdx] = 0.0;
        values[temperatureIdx] = 283.0 + globalPos[1]*0.03;
    }
};
} //end namespace

#endif
