// $Id$
#ifndef DUNE_FETOTALVELOCITY2P_HH
#define DUNE_FETOTALVELOCITY2P_HH
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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

#include "dumux/functions/CRfunction.hh"

namespace Dune
{

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
void FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, Communication>::calculateVelocity(const Scalar t=0) const
{
	typedef typename GridView::Grid Grid;
	enum {dim=Grid::dimension};
	typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
	typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
	typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;

	ElementMapper elementMapper(this->gridView);

	ElementIterator endEIt = this->gridView.template end<0>();
	for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != endEIt; ++eIt)
	{
		const Element& element = *eIt;

        // element geometry
        const Geometry& geometry = element.geometry();

		Dune::GeometryType geomType = geometry.type();

		const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
		Dune::FieldVector<Scalar,dim> global = geometry.global(local);

		const typename CRShapeFunctionSetContainer<double,double,dim>::value_type&
		sfs = CRShapeFunctions<double,double,dim>::general(geomType,1);

		int eIdx = elementMapper.map(element);

        // get the absolute permeability
		FieldMatrix<double,dim,dim> K = this->diffProblem.soil().K(global, element, local);

		IntersectionIterator endis = this->gridView.iend(element);
		for (IntersectionIterator is = this->gridView.ibegin(element); is!=endis; ++is)
		{
			// get geometry type of face
			GeometryType gtf = is->geometryInInside().type();

			// local number of facet
			int i = is->indexInInside();

			const FieldVector<double,dim>& faceLocal = sfs[i].position();

			FieldVector<double,dim> gradP(0);
			for (int comp = 0; comp < dim; comp++)
			{
				FieldVector<int,dim> order(0);
				order[comp] = 1;
				gradP[comp] = pressP1.derivativelocal(comp, order, element, faceLocal);
			}

			// get the negative exact velocity
			FieldVector<double,dim> KGradP(0);
			K.umv(gradP, KGradP);

			KGradP *= -1.0;

			this->diffProblem.variables().velocity()[eIdx][i] = KGradP;
		}
	}

	//DUNE_THROW(Dune::NotImplemented, "velocities only implemented for finite volume and mimetic finite differences discretisations");
    return;
}

}
#endif
