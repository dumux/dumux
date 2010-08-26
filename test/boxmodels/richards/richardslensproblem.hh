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
#ifndef DUMUX_RICHARDS_LENSPROBLEM_HH
#define DUMUX_RICHARDS_LENSPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/boxmodels/richards/richardsmodel.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "richardslensspatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class RichardsLensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(RichardsLensProblem, INHERITS_FROM(BoxRichards));

// Set the grid type
// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(RichardsLensProblem, Grid, Dune::UGGrid<2>);
#else
SET_PROP(RichardsLensProblem, Grid) { typedef Dune::SGrid<2, 2> type; };
//SET_TYPE_PROP(RichardsLensProblem, Grid, Dune::YaspGrid<2>);
#endif

SET_PROP(RichardsLensProblem, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for cubes
//    typedef Dune::PDELab::P1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for simplices
};

// Set the problem property
SET_PROP(RichardsLensProblem, Problem)
{
    typedef Dumux::RichardsLensProblem<TTAG(RichardsLensProblem)> type;
};

// Set the wetting phase
SET_PROP(RichardsLensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(RichardsLensProblem, SpatialParameters)
{
    typedef Dumux::RichardsLensSpatialParameters<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(RichardsLensProblem, EnableGravity, true);

// Write the intermediate results of the newton method?
SET_BOOL_PROP(RichardsLensProblem, NewtonWriteConvergence, false);
}

/*!
 * \ingroup RichardsBoxProblems
 * \brief Base class for defining an instance of a Richard`s problem, where water is infiltrating into an initially unsaturated zone.
 *
 * The domain is box shaped. Left and right boundaries are Dirichlet boundaries with fixed water pressure (fixed Saturation Sw = 0),
 * bottom boundary is closed (Neumann 0 boundary), the top boundary (Neumann 0 boundary) is also closed except for infiltration section,
 * where water is infiltrating into an initially unsaturated zone. Linear capillary pressure-saturation relationship is used.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./lens_richards ./grids/lens_richards_2d.dgf 10000 1</tt>
 *
 * where start simulation time = 1 second, end simulation time = 10000 seconds
 */
template <class TypeTag = TTAG(RichardsLensProblem) >
class RichardsLensProblem : public RichardsBoxProblem<TypeTag>
{
    typedef RichardsLensProblem<TypeTag> ThisType;
    typedef RichardsBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        // copy some indices for convenience
        pwIdx = Indices::pwIdx,
        contiEqIdx = Indices::contiEqIdx,

        // Grid and world dimension
        dimWorld = GridView::dimensionworld,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    RichardsLensProblem(TimeManager &timeManager,
                        const GridView &gridView,
                        const GlobalPosition &lensLowerLeft,
                        const GlobalPosition &lensUpperRight)
        : ParentType(timeManager, gridView)
    {
        lensLowerLeft_=lensLowerLeft;
        lensUpperRight_=lensUpperRight;
        this->spatialParameters().setLensCoords(lensLowerLeft_, lensUpperRight_);
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
    { return "richardslens"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
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
        return 1e5; // reference non-wetting phase pressure [Pa] used for viscosity and density calculations
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
    void boundaryTypes(BoundaryTypes &values,
                       const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       const Intersection &is,
                       int scvIdx,
                       int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        if (onLeftBoundary_(globalPos) || 
            onRightBoundary_(globalPos))
        {
            values.setAllDirichlet();
        }
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values,
                   const Element &element,
                   const FVElementGeometry &fvElemGeom,
                   const Intersection &is,
                   int scvIdx,
                   int boundaryFaceIdx) const
    {
        // use initial values as boundary conditions
        initial(values, element, fvElemGeom, scvIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = 0.0;
        if (onInlet_(globalPos)) {
            // inflow of water
            values[contiEqIdx] = -0.04; // kg / (m * s)
        }
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
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        //const GlobalPosition &pos = element.geometry().corner(scvIdx);
        
        Scalar Sw = 0.0;
        Scalar pc =
            MaterialLaw::pC(this->spatialParameters().materialLawParams(element, 
                                                                        fvElemGeom,
                                                                        scvIdx),
                            Sw);
        values[pwIdx] = pNreference() - pc;
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
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    static const Scalar eps_ = 3e-6;
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
};
} //end namespace

#endif
