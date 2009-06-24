// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Melanie Darcis                                    *
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
#ifndef DUNE_INJECTIONPROBLEM2PNI_HH
#define DUNE_INJECTIONPROBLEM2PNI_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/matrixproperties.hh>
#include <dumux/material/fluids/air.hh>
#include <dumux/material/fluids/water.hh>
#include <dumux/boxmodels/2pni/2pniboxmodel.hh>

namespace Dune {

template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(BoxTwoPNI));

// Set the grid type
SET_PROP(InjectionProblem2PNI, Grid)
{
#if ENABLE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::SGrid<2, 2> type;
    //typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(InjectionProblem2PNI, Problem)
{
    typedef Dune::InjectionProblem2PNI<TTAG(InjectionProblem2PNI)> type;
};

// Set the wetting phase
SET_TYPE_PROP(InjectionProblem2PNI, WettingPhase, Dune::Water);

// Set the non-wetting phase
SET_TYPE_PROP(InjectionProblem2PNI, NonwettingPhase, Dune::Air);

// Set the soil properties
SET_PROP(InjectionProblem2PNI, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dune::HomogeneousSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem2PNI, EnableGravity, true);
}

/*!
 * \ingroup TwoPNIBoxProblems
 * \brief Nonisothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area produced by a heat source.
 *
 * The domain is sized 40 m times 40 m. The rectangular heat source starts at (20 m, 5 m)
 * and ends at (30 m, 35 m)
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top and on the bottom of the domain, while dirichlet conditions
 * apply on the left and the right boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the bottom boundary from 15 m to 25 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * This problem uses the \ref TwoPNIBoxModel.
 *
 * This problem should typically be simulated for 300000 s is reached.
 * A good choice for the initial time step size is 1000 s.
 */
template<class TypeTag = TTAG( InjectionProblem2PNI)>
class InjectionProblem2PNI
  : public TwoPNIBoxProblem<TypeTag,
                            InjectionProblem2PNI<TypeTag> >
{
    typedef InjectionProblem2PNI<TypeTag>                     ThisType;
    typedef TwoPNIBoxProblem<TypeTag, ThisType>               ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPNIIndices)) Indices;
	enum {
		numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
		pressureIdx = Indices::pressureIdx,
		saturationIdx = Indices::saturationIdx,
		temperatureIdx = Indices::temperatureIdx,

		// Grid and world dimension
		dim = GridView::dimension,
		dimWorld = GridView::dimensionworld,
	};

	typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
	typedef typename SolutionTypes::PrimaryVarVector PrimaryVarVector;
	typedef typename SolutionTypes::BoundaryTypeVector BoundaryTypeVector;

	typedef typename GridView::template Codim<0>::Entity Element;
	typedef typename GridView::template Codim<dim>::Entity Vertex;
	typedef typename GridView::IntersectionIterator IntersectionIterator;

	typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

	typedef Dune::FieldVector<Scalar, dim> LocalPosition;
	typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
	InjectionProblem2PNI(const GridView &grid)
        : ParentType(grid)
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
    { return "injection2pni"; }

    // \}

	/*!
	 * \name Boundary conditions
	 */
	// \{

	/*!
	 * \brief Specifies which kind of boundary condition should be
	 *        used for which equation on a given boundary segment.
	 */
	void boundaryTypes(BoundaryTypeVector &values,
			const Element &element,
			const FVElementGeometry &fvElemGeom,
			const IntersectionIterator &isIt,
			int scvIdx,
			int boundaryFaceIdx) const
	{
		const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        if(globalPos[0] > 40 - eps_ || globalPos[0] < eps_)
            values = BoundaryConditions::dirichlet;
        else
            values = BoundaryConditions::neumann;

#if !ISOTHERMAL
        values[temperatureIdx] = BoundaryConditions::dirichlet;
#endif
	}

	/*!
	 * \brief Evaluate the boundary conditions for a dirichlet
	 *        boundary segment.
	 *
	 * For this method, the \a values parameter stores primary variables.
	 */
	void dirichlet(PrimaryVarVector &values,
			const Element &element,
			const FVElementGeometry &fvElemGeom,
			const IntersectionIterator &isIt,
			int scvIdx,
			int boundaryFaceIdx) const
	{
		const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

	         Scalar densityW = 1000.0;
	         values[pressureIdx] = 1e5 + (depthBOR_ - globalPos[1])*densityW*9.81;
	         values[saturationIdx] = 0.0;
	 #if !ISOTHERMAL
	         values[temperatureIdx] = 283.0 + (depthBOR_ - globalPos[1])*0.03;
	 #endif
	}

	/*!
	 * \brief Evaluate the boundary conditions for a neumann
	 *        boundary segment.
	 *
	 * For this method, the \a values parameter stores the mass flux
	 * in normal direction of each phase. Negative values mean influx.
	 */
	void neumann(PrimaryVarVector &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int scvIdx,
                 int boundaryFaceIdx) const
	{
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        values = 0;
        // negative values for injection
        if (globalPos[0] > 15 && globalPos[0] < 25 &&
            globalPos[1] < eps_)
        {
            values[saturationIdx] = -1e-3;
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
	void source(PrimaryVarVector &values,
			const Element &element,
			const FVElementGeometry &,
			int subControlVolumeIdx) const
	{
	       const GlobalPosition &globalPos = element.geometry().corner(subControlVolumeIdx);
	       values = Scalar(0.0);
	       if (globalPos[0] > 20 && globalPos[0] < 30 && globalPos[1] > 5 && globalPos[1] < 35)
	    	  values[temperatureIdx] = 100.;
	}

	/*!
	 * \brief Evaluate the initial value for a control volume.
	 *
	 * For this method, the \a values parameter stores primary
	 * variables.
	 */
	void initial(PrimaryVarVector &values,
			const Element &element,
			const FVElementGeometry &fvElemGeom,
			int scvIdx) const
	{
		const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

	         Scalar densityW = 1000.0;
	         values[pressureIdx] = 1e5 + (depthBOR_ - globalPos[1])*densityW*9.81;
	         values[saturationIdx] = 0.0;
	 #if !ISOTHERMAL
	         values[temperatureIdx] = 283.0 + (depthBOR_ - globalPos[1])*0.03;
		     if (globalPos[0] > 20 && globalPos[0] < 30 && globalPos[1] > 5 && globalPos[1] < 35)
		    	values[temperatureIdx] = 380;
	 #endif
	}
	// \}

private:

	static const Scalar eps_ = 1e-6;
	static const Scalar depthBOR_ = 1000.0;
};
} //end namespace

#endif
