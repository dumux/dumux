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
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 */

#ifndef DUMUX_INJECTION_PROBLEM_2PNI_HH
#define DUMUX_INJECTION_PROBLEM_2PNI_HH

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>

#include "injection2pspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(TwoPNI, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionCCProblem2PNI, INHERITS_FROM(CCModel, InjectionProblem2PNI));

// Set the grid type
SET_TYPE_PROP(InjectionProblem2PNI, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InjectionProblem2PNI, Problem, InjectionProblem2PNI<TypeTag>);

// Use the same fluid system as the 2p2c injection problem
SET_TYPE_PROP(InjectionProblem2PNI, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 *
 * The domain is sized 60 m times 40 m. The rectangular area with the increased temperature (380 K)
 * starts at (20 m, 5 m) and ends at (30 m, 35 m)
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top, on the bottom and on the right of the domain, while dirichlet conditions
 * apply on the left boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the right boundary from 7 m to 15 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * This problem uses the \ref TwoPModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./exercise1_2pni -ParameterFile exercise1.input</tt>
 */
template<class TypeTag>
class InjectionProblem2PNI : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        contiNEqIdx = Indices::contiNEqIdx,

    /*
     * dumux-course-task:
     * get the temperatureIdx as the index for the primary variable temperature and the
     * energyIdx as the index for the energy conservation equation for your convinience.
     */

        // world dimension
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem2PNI(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        maxDepth_  = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.MaxDepth);

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                /*tempMax=*/423.15,
                /*numTemp=*/50,
                /*pMin=*/0.0,
                /*pMax=*/30e6,
                /*numP=*/300);
    }

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
    { return "injection-2pni";}

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
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

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

         /*!
          * dumux-course-task:
          * set dirichlet conditions for the temperature index everywhere
          */
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        values[saturationIdx] = 0.0;
       /*!
        * dumux-course-task:
        * set a temperature gradient of 0.03 K per m beginning at 283 K here.
        * Hint: you can use maxDepth_ and the globalPos similar to the pressure gradient
        */
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     *
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        values = 0;

        GlobalPosition globalPos;
        if (isBox)
            globalPos = element.geometry().corner(scvIdx);
        else
            globalPos = intersection.geometry().center();

        if (globalPos[1] < 15 + eps_ && globalPos[1] > 7 - eps_) {
            // inject air. negative values mean injection
            values[contiNEqIdx] = -1e-3; // kg/(s*m^2)
        }
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
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        values[saturationIdx] = 0.0;

        /*!
        * dumux-course-task:
        * set a temperature gradient of 0.03 K per m beginning at 283 K here.
        * Hint: you can use maxDepth_ and the globalPos similar to the pressure gradient
        * use globalPos[0] and globalpos[1] to implement the high temperature lens with 380 K
        */
    }


private:
    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};
} //end namespace

#endif
