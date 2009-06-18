// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_RICHARDSPROBLEM_HH
#define DUNE_RICHARDSPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/air.hh>

#include <dumux/boxmodels/richards/richardsboxmodel.hh>

#include "richardssoil.hh"

namespace Dune
{

template <class TypeTag>
class RichardsTestProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(RichardsTestProblem, INHERITS_FROM(BoxRichards));

// Set the grid type
// Set the grid type
#if ENABLE_UG
SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::UGGrid<2>);
#else
SET_PROP(RichardsTestProblem, Grid) { typedef Dune::SGrid<2, 2> type; };
//SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_PROP(RichardsTestProblem, Problem)
{
    typedef Dune::RichardsTestProblem<TTAG(RichardsTestProblem)> type;
};

// Set the wetting phase
SET_TYPE_PROP(RichardsTestProblem, WettingPhase, Dune::Water);

// Set the non-wetting phase
SET_TYPE_PROP(RichardsTestProblem, NonwettingPhase, Dune::Air);

// Set the soil properties
SET_PROP(RichardsTestProblem, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    
public:
    typedef Dune::RichardsSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(RichardsTestProblem, EnableGravity, true);
}

/*!
 * \ingroup RichardsBoxProblems
 * \brief  Base class for defining an instance of a Richard`s problem, where water is inflitrating into an initially unsaturated zone.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary) except the top boundary(Dirichlet),
 * where water is inflitrating into an initially unsaturated zone. Linear capillary pressure-saturation relationship is used.
 *
 * To run the simulation execute the following line in shell:
 * ./new_test_richards ./grids/richards_2d.dgf 10000 1
 * where start simulation time = 1 second, end simulation time = 10000 seconds
 * The same file can be also used for 3d simulation but you need to change line
 * "typedef Dune::UGGrid<2> type;" with "typedef Dune::UGGrid<3> type;" and use richards_3d.dgf grid
 */
template <class TypeTag = TTAG(RichardsTestProblem) >
class RichardsTestProblem : public RichardsBoxProblem<TypeTag, 
                                                      RichardsTestProblem<TypeTag> >
{
    typedef RichardsTestProblem<TypeTag>   ThisType;
    typedef RichardsBoxProblem<TypeTag, ThisType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // copy some indices for convenience
        pW = Indices::pW,

        // Grid and world dimension
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;

    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef typename GridView::template Codim<dim>::Entity      Vertex;
    typedef typename GridView::IntersectionIterator             IntersectionIterator;
  
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

public:
    RichardsTestProblem(const GridView &gridView)
        : ParentType(gridView)
    {
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
    const char *name() const
    { return "richardstest"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 283.15; // -> 10°C
    };

    /*!
     * \brief Returns the reference pressure of the nonwetting phase.
     *
     * This problem assumes a pressure of 1 bar.
     */
    Scalar pNreference() const
    {
        return 1.0e+5; // reference non-wetting phase pressure [Pa] used for viscosity and density calculations
    };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*! 
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;

        if (globalPos[dim-1] > this->bboxMax()[dim-1] - eps)
            // dirichlet on the lower boundary
            values = BoundaryConditions::dirichlet;
        else
            values = BoundaryConditions::neumann;
    }

    /*! 
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVarVector           &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
    	const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        
        if (globalPos[dim-1] > this->bboxMax()[dim-1] - eps)
        	values[pW] = 1.0e+5*0.99;
        else 
            values[pW] = 0.0;
    }

    /*! 
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVarVector           &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        // const GlobalPosition &globalPos = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
    	values[pW] = 0;
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*! 
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate
     * wetting phase mass is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*! 
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVarVector        &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        // const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values[pW] = 1.0e+5 - 5.0e+4;
    }
    // \}

private:
    static const Scalar eps = 1e-6;
};
} //end namespace

#endif
