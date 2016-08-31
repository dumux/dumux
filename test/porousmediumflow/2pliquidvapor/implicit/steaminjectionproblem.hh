// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later vesion.                                      *
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
 * \brief Non-isothermal steam injection problem for the validation of the 2p1cni model against the data provided by Ochs, 2006.
 *
 */
#ifndef DUMUX_STEAM_INJECTIONPROBLEM_HH
#define DUMUX_STEAM_INJECTIONPROBLEM_HH

#include <dumux/porousmediumflow/2pliquidvapor/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "steaminjectionspatialparams.hh"
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>

namespace Dumux
{
template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem, INHERITS_FROM(TwoPOneCNI, InjectionProblemSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem, INHERITS_FROM(BoxModel, InjectionProblem));
NEW_TYPE_TAG(InjectionCCProblem, INHERITS_FROM(CCModel, InjectionProblem));

SET_TYPE_PROP(InjectionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_PROP(InjectionProblem, Problem)
{
    typedef Dumux::InjectionProblem<TypeTag> type;
};

// Set the fluid system
SET_PROP(InjectionProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> >  H2OType ;
public:
    typedef Dumux::FluidSystems::TwoPLiquidVaporFluidsystem<Scalar, H2OType > type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem, ProblemEnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(InjectionProblem, ImplicitNumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(InjectionProblem, NewtonWriteConvergence, false);

//Define whether spurious cold-water flow into the steam is blocked
SET_BOOL_PROP(InjectionProblem, UseBlockingOfSpuriousFlow, true);
}


/*!
 * \ingroup ThreePTwoCNIBoxModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal problem where ...
 *
 * This problem uses the \ref ThreePTwoCNIModel.
 *
 *  */
template <class TypeTag >
class InjectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,

        energyEqIdx = Indices::energyEqIdx,

        // phase and component indices
        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,


        // Phase State
        wPhaseOnly = Indices::wPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        twoPhases = Indices::twoPhases,

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
    typedef Dune::FieldVector<Scalar,dim> Vector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;


    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */

    InjectionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-6)
    {

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        FluidSystem::init();
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
    { return name_; }


    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
     {
          values = 0.0;
     }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
   void boundaryTypesAtPos(BoundaryTypes &bcTypes,
            const GlobalPosition &globalPos) const
    {
        if(globalPos[1] > this->bBoxMax()[1] - eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
           bcTypes.setAllDirichlet();
        else
           bcTypes.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
       initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeomtry The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */

      void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;

        if (globalPos[0] < eps_)
        {
            if(globalPos[1] > 2.0 - eps_ && globalPos[1] < 3.0 + eps_)
            {
            Scalar massRate = 1e-1;
            values[Indices::conti0EqIdx] = -massRate;
            values[Indices::energyEqIdx] = -massRate * 2690e3;
            }
        }
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
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vert The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vert,
                             int &globalIdx,
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
        values[pressureIdx] = 101300.0 + (this->bBoxMax()[1] - globalPos[1])*densityW*9.81; // hydrostatic pressure
        values[switch1Idx] = 283.13;
    }

    Scalar eps_;
    std::string name_;
};
} //end namespace

#endif
