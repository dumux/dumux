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

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dumux/io/simplexgridcreator.hh>
#include <dumux/io/cubegridcreator.hh>

#include <dumux/implicit/2p/2pmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>

// use the same spatial parameters as the injection problem of the
// 2p2c test program
#include "../2p2c/injectionspatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux {

template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
#if !ISOTHERMAL
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(TwoPNI, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem2PNI, INHERITS_FROM(BoxModel, InjectionProblem2PNI));
NEW_TYPE_TAG(InjectionCCProblem2PNI, INHERITS_FROM(CCModel, InjectionProblem2PNI));
#else
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(TwoP, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem2PNI, INHERITS_FROM(BoxModel, InjectionProblem2PNI));
NEW_TYPE_TAG(InjectionCCProblem2PNI, INHERITS_FROM(CCModel, InjectionProblem2PNI));
#endif

// set the GridCreator property
#if HAVE_UG
SET_TYPE_PROP(InjectionProblem2PNI, GridCreator, Dumux::SimplexGridCreator<TypeTag>);
#else
SET_TYPE_PROP(InjectionProblem2PNI, GridCreator, Dumux::CubeGridCreator<TypeTag>);
#endif

// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(InjectionProblem2PNI, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(InjectionProblem2PNI, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(InjectionProblem2PNI, Problem, Dumux::InjectionProblem2PNI<TypeTag>);

#if 1
// Use the same fluid system as the 2p2c injection problem
SET_TYPE_PROP(InjectionProblem2PNI, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);
#else
// Set the wetting phase
SET_TYPE_PROP(InjectionProblem2PNI,
              WettingPhase,
              Dumux::LiquidPhase<GET_PROP_TYPE(TypeTag, Scalar) Scalar,
                                 Dumux::SimpleH2O<GET_PROP_TYPE(TypeTag, Scalar) Scalar> >);

// Set the non-wetting phase
SET_TYPE_PROP(InjectionProblem2PNI,
              NonwettingPhase,
              Dumux::GasPhase<GET_PROP_TYPE(TypeTag, Scalar) Scalar,
                                 Dumux::N2<GET_PROP_TYPE(TypeTag, Scalar) Scalar> >);
#endif
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
 * Gas is injected at the right boundary from 5 m to 15 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * This problem uses the \ref TwoPModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pni -parameterFile test_box2pni.input</tt> or
 * <tt>./test_cc2pni -parameterFile test_cc2pni.input</tt>
 */
template<class TypeTag>
class InjectionProblem2PNI : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

#if ISOTHERMAL
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
#else
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
#endif
    enum {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        contiNEqIdx = Indices::contiNEqIdx,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

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
        maxDepth_ = 2700.0; // [m]
        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                /*tempMax=*/423.15,
                /*numTemp=*/50,
                /*pMin=*/0.0,
                /*pMax=*/30e6,
                /*numP=*/300);
        
        name_               = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
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
    { return name_; }

#if ISOTHERMAL
    /*!
     * \brief Returns the temperature \f$ K \f$
     */
    Scalar temperature() const
    {
        return 273.15 + 30; // [K]
    }
#endif


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

#if !ISOTHERMAL
        // set a dirichlet value for the temperature, use the energy
        // equation to set the value
        values.setDirichlet(temperatureIdx, energyEqIdx);
#endif
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
#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;
#endif
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
        
        if (globalPos[1] < 15 && globalPos[1] > 7) {
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

#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;
        if (globalPos[0] > 20 && globalPos[0] < 30 && globalPos[1] > 5 && globalPos[1] < 35)
            values[temperatureIdx] = 380;
#endif // !ISOTHERMAL
    }
    // \}

private:
    Scalar maxDepth_;
    Scalar eps_;
    std::string name_;
};
} //end namespace

#endif
