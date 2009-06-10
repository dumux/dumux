//$Id$
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

#include <dumux/material/fluids/brine_co2.hh>
#include <dumux/boxmodels/2pni/2pniboxmodel.hh>

#include "injectionsoil2pni.hh"


namespace Dune {

template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(BoxTwoPNI));

// Set the grid type
SET_TYPE_PROP(InjectionProblem2PNI, Grid, Dune::UGGrid<2>);

// Set the problem property
SET_PROP(InjectionProblem2PNI, Problem)
{
    typedef Dune::InjectionProblem2PNI<TTAG(InjectionProblem2PNI)> type;
};

// Set the wetting phase
SET_TYPE_PROP(InjectionProblem2PNI, WettingPhase, Dune::Liq_BrineCO2);

// Set the non-wetting phase
SET_TYPE_PROP(InjectionProblem2PNI, NonwettingPhase, Dune::Gas_BrineCO2);

// Set the soil properties
SET_PROP(InjectionProblem2PNI, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    
public:
    typedef Dune::InjectionSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem2PNI, EnableGravity, true);
}

/*!
 * \ingroup TwoPNIBoxProblems
 * \brief Gas injection problem where a gas (e.g. \f$ CO_2 \f$ is injected into a fully
 *        water saturated medium.
 *
 * The domain is sized 100m times 50m and features a rectangular lens
 * with low permeablility which spans from (0 m , 30 m) to (50 m, 35 m)
 * and is surrounded by a medium with higher permability.
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the left, the top and the bottom of the domain, while dirichlet conditions
 * apply on the right boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the left boundary from 0 m to 20 m at a rate of
 * 0.04 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature are applied.
 *
 * This problem uses the \ref TwoPNIBoxModel.
 *
 * This problem should typically simulated until \f$ t_{\text{end}} =
 * 50\,000\;s \f$ is reached. A good choice for the initial time step
 * size is \f$ t_{\text{inital}} = 1\,000\;s \f$
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
	InjectionProblem2PNI(const GridView &grid,
                         const GlobalPosition &lensLowerLeft,
                         const GlobalPosition &lensUpperRight)
        : ParentType(grid)
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

		values = BoundaryConditions::neumann;

		if (onRightBoundary_(globalPos))
            values = BoundaryConditions::dirichlet;

		values[temperatureIdx] = BoundaryConditions::dirichlet;
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

		Scalar densityB = 1046;
		Scalar pRef = 101300.;
		Scalar TRef = 283.15;

		values[pressureIdx] = pRef - (depthBOR_ - globalPos[dim - 1]) * densityB*this->gravity()[dim-1];
		values[saturationIdx] = 0.0;
		values[temperatureIdx] = TRef + (depthBOR_ - globalPos[dim - 1]) * 0.03;
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
		if (onInlet_(globalPos))
            values[Indices::phase2Mass(Indices::nPhase)] = -0.0004; // kg / (m * s)
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
		values = Scalar(0.0);
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

		Scalar densityB = 1046;
		Scalar pRef = 101300.;
		Scalar TRef = 283.15;

		values[pressureIdx] = pRef - (depthBOR_ - globalPos[dim-1]) * densityB*this->gravity()[dim-1];
		values[saturationIdx] = 0.0;
		values[temperatureIdx] = TRef + (depthBOR_ - globalPos[dim-1]) * 0.03;
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

	bool onInlet_(const GlobalPosition &globalPos) const
	{
		Scalar height = this->bboxMax()[1] - this->bboxMin()[1];
		Scalar lambda = globalPos[1]/height;
		return onLeftBoundary_(globalPos) && 0.0 < lambda && lambda < 1.0/5.0;
	}

	static const Scalar eps_ = 1e-6;
	static const Scalar depthBOR_ = 2500.0;
};
} //end namespace

#endif
