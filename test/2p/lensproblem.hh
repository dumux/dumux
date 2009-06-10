/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUNE_LENSPROBLEM_HH
#define DUNE_LENSPROBLEM_HH

#define USE_UG 1

#if USE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif

#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/dnapl.hh>

#include <dumux/boxmodels/2p/2pboxmodel.hh>

#include "lenssoil.hh"

namespace Dune
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(BoxTwoP));

// Set the grid type
SET_PROP(LensProblem, Grid)
{
#if USE_UG
    typedef Dune::UGGrid<2> type;
#else // USE_UG
    typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(LensProblem, Problem)
{
    typedef Dune::LensProblem<TTAG(LensProblem)> type;
};

// Set the wetting phase
SET_TYPE_PROP(LensProblem, WettingPhase, Dune::Water);

// Set the non-wetting phase
SET_TYPE_PROP(LensProblem, NonwettingPhase, Dune::DNAPL);

// Set the soil properties
SET_PROP(LensProblem, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    
public:
    typedef Dune::LensSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(LensProblem, EnableGravity, true);
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief Soil decontamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability.
 * 
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically simulated until \f$t_{\text{end}} =
 * 50\,000\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\,000\;s\f$
 */
template <class TypeTag = TTAG(LensProblem) >
class LensProblem  : public TwoPBoxProblem<TypeTag, 
                                           LensProblem<TypeTag> >
{
    typedef LensProblem<TypeTag>   ThisType;
    typedef TwoPBoxProblem<TypeTag, ThisType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // copy some indices for convenience
        pressureIdx   = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pW = Indices::pW,
        sN = Indices::sN,

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
    LensProblem(const GridView &gridView,
                const GlobalPosition &lensLowerLeft,
                const GlobalPosition &lensUpperRight)
        : ParentType(gridView)
    {
        this->soil().setLensCoords(lensLowerLeft, lensUpperRight);
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
    { return "lens"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 283.15; // -> 10°C
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
            = element.geometry().corner(scvIdx);

        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            values[pW] = BoundaryConditions::dirichlet;
            values[sN] = BoundaryConditions::dirichlet;
        }
        else {
            values[pW] = BoundaryConditions::neumann;
            values[sN] = BoundaryConditions::neumann;
        }
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        
        Scalar densityW = this->wettingPhase().density();
        
        if (onLeftBoundary_(globalPos))
        {
            Scalar height = this->bboxMax()[1] - this->bboxMin()[1];
            Scalar depth = this->bboxMax()[1] - globalPos[1];
            Scalar alpha = (1 + 0.5/height);

            // hydrostatic pressure scaled by alpha
            values[pW] = - alpha*densityW*this->gravity()[1]*depth;
            values[sN] = 0.0;
        }
        else if (onRightBoundary_(globalPos))
        {
            Scalar depth = this->bboxMax()[1] - globalPos[1];

            // hydrostatic pressure
            values[pW] = -densityW*this->gravity()[1]*depth;
            values[sN] = 0.0;
        }
        else
            values = 0.0;
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = 0.0;
        if (onInlet_(globalPos)) {
            values[Indices::phase2Mass(Indices::nPhase)] = -0.04; // kg / (m * s)
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*! 
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
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
        // no DNAPL, some random pressure
        values[pW] = 0.0;
        values[sN] = 0.0;
    }
    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bboxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bboxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bboxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bboxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        Scalar lambda = (this->bboxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda  && lambda < 2.0/3.0;
    }
    
    static const Scalar eps_ = 3e-6;
};
} //end namespace

#endif
