// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
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

#ifndef DUNE_MIMETICGROUNDWATER_HH
#define DUNE_MIMETICGROUNDWATER_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>

#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/groundwater/groundwater.hh>
#include"dumux/shapefunctions/CRshapefunctions.hh"

/**
* @file
* @brief  compute local stiffness matrix for conforming finite elements for diffusion equation
* @author Peter Bastian
*/

namespace Dune
{
/** @addtogroup DISC_Disc
*
* @{
*/
/**
* @brief compute local stiffness matrix for conforming finite elements for diffusion equation
*
*/

//! A class for computing local stiffness matrices
/*! A class for computing local stiffness matrix for the
diffusion equation

div j = q; j = -K grad u; in Omega

u = g on Gamma1; j*n = J on Gamma2.

Uses conforming finite elements with the CR shape functions.
It should work for all dimensions and element types.
All the numbering is with respect to the reference element and the
CR shape functions

Template parameters are:

- Grid  a DUNE grid type
- RT    type used for return values
*/
template<class GridView, class Scalar, class VC, class Problem>
class MimeticGroundwaterEquationLocalStiffness
:
public LocalStiffness<GridView, Scalar, 1>
{
	// grid types
	enum
	{
		dim = GridView::dimension
	};
	typedef typename GridView::Grid Grid;
	typedef typename GridView::Traits::template Codim<0>::Entity Element;

public:
	// define the number of components of your system, this is used outside
	// to allocate the correct size of (dense) blocks with a FieldMatrix
	enum
	{   m=1};
	enum
	{   size=CRShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};

	//! Constructor
	MimeticGroundwaterEquationLocalStiffness (Problem& problem,
			bool levelBoundaryAsDirichlet, const GridView& gridView,
			bool procBoundaryAsDirichlet=true)
	: problem_(problem),levelBoundaryAsDirichlet_(levelBoundaryAsDirichlet),
	procBoundaryAsDirichlet_(procBoundaryAsDirichlet), gridView_(gridView)
	{}

	//! assemble local stiffness matrix for given element and order
	/*! On exit the following things have been done:
	- The stiffness matrix for the given entity and polynomial degree has been assembled and is
	accessible with the mat() method.
	- The boundary conditions have been evaluated and are accessible with the bc() method
	- The right hand side has been assembled. It contains either the value of the essential boundary
	condition or the assembled source term and neumann boundary condition. It is accessible via the rhs() method.
	@param[in]  e    a codim 0 entity reference
	@param[in]  k    order of CR basis
	*/
	void assemble (const Element& element, int k=1)
	{
		// extract some important parameters
		Dune::GeometryType gt = element.geometry().type();
		const typename Dune::CRShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
		sfs=Dune::CRShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
		setcurrentsize(sfs.size());

		// clear assemble data
		for (int i=0; i<sfs.size(); i++)
		{
			this->b[i] = 0;
			this->bctype[i][0] = BoundaryConditions::neumann;
			for (int j=0; j<sfs.size(); j++)
				this->A[i][j] = 0;
		}

		assembleV(element,k);
		assembleBC(element,k);
	}

	// TODO/FIXME: this is only valid for linear problems where
	// the local stiffness matrix is independend of the current
	// solution. We need to implement this properly, but this
	// should at least make the thing compile...
	typedef Dune::FieldVector<Scalar, m> VBlockType;
	void assemble(const Element &cell, const Dune::BlockVector<VBlockType>& localSolution, int orderOfShapeFns = 1)
	{
		assemble(cell, orderOfShapeFns);
	}

	//! assemble only boundary conditions for given element
	/*! On exit the following things have been done:
	- The boundary conditions have been evaluated and are accessible with the bc() method
	- The right hand side contains either the value of the essential boundary
	condition or the assembled neumann boundary condition. It is accessible via the rhs() method.
	@param[in]  element    a codim 0 entity reference
	@param[in]  k    order of CR basis
	*/
	void assembleBoundaryCondition (const Element& element, int k=1)
	{
		// extract some important parameters
		Dune::GeometryType gt = element.geometry().type();
		const typename Dune::CRShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
		sfs=Dune::CRShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
		setcurrentsize(sfs.size());

		// clear assemble data
		for (int i=0; i<sfs.size(); i++)
		{
			this->b[i] = 0;
			this->bctype[i][0] = BoundaryConditions::neumann;
		}

		this->template assembleBC(element,k);
	}

	void assembleElementMatrices(const Element& element, Dune::FieldVector<Scalar,2*dim>& faceVol,
			Dune::FieldMatrix<Scalar,2*dim,2*dim>& W, Dune::FieldVector<Scalar,2*dim>& c,
			Dune::FieldMatrix<Scalar,2*dim,2*dim>& Pi, Scalar& dinv, Dune::FieldVector<Scalar,2*dim>& F, Scalar& qmean)
	{
		// extract some important parameters
		Dune::GeometryType gt = element.geometry().type();
		const typename Dune::CRShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
		sfs=Dune::CRShapeFunctions<Scalar,Scalar,dim>::general(gt,1);
		setcurrentsize(sfs.size());

		// cell center in reference element
		const Dune::FieldVector<Scalar,dim>& centerLocal = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0,0);

		// get global coordinate of cell center
		Dune::FieldVector<Scalar,dim> centerGlobal = element.geometry().global(centerLocal);

		int globalIdx = problem_.variables().indexDiffusion(element);

		// eval diffusion tensor, ASSUMING to be constant over each cell
		Dune::FieldMatrix<Scalar,dim,dim> K(0);
		K = problem_.soil().K(centerGlobal,element,centerLocal);

		K *= (problem_.variables().mobilityWetting()[globalIdx] + problem_.variables().mobilityNonWetting()[globalIdx]);

		// cell volume
		Scalar volume = element.geometry().volume();

		// build the matrices R and ~N
		Dune::FieldMatrix<Scalar,2*dim,dim> R(0), N(0);

		//       std::cout << "element " << elemId << ": center " << centerGlobal << std::endl;;

		typedef typename GridView::IntersectionIterator IntersectionIterator;

		IntersectionIterator endit = gridView_.template iend(element);
		for (IntersectionIterator it = gridView_.template ibegin(element); it!=endit; ++it)
		{
			// get geometry type of face
			Dune::GeometryType gtf = it->geometryInInside().type();

			// local number of facet
			int i = it->indexInInside();

			const Dune::FieldVector<Scalar,dim>& faceLocal = sfs[i].position();
			Dune::FieldVector<Scalar,dim> faceGlobal = element.geometry().global(faceLocal);
			faceVol[i] = it->geometry().volume();

			//           std::cout << "  face " << i << ": local = " << faceLocal << ", global = " << faceGlobal << std::endl;
			//           std::cout << "    boundary = " << it.boundary() << ", neighbor = " << it->neighbor() << std::endl;
			//", outside elemId = " << elementmapper.map(*(it->outside())) << std::endl;

			// center in face's reference element
			const Dune::FieldVector<Scalar,dim-1>&
			faceLocalNm1 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);

			// get normal vector
			Dune::FieldVector<Scalar,dim> unitOuterNormal = it->unitOuterNormal(faceLocalNm1);

			N[i] = unitOuterNormal;

			for (int k = 0; k < dim; k++)
				// move origin to the center of gravity
				R[i][k] = faceVol[i]*(faceGlobal[k] - centerGlobal[k]);
		}

		//      std::cout << "N =\dim" << N;
		//      std::cout << "R =\dim" << R;

		// proceed along the lines of Algorithm 1 from
		// Brezzi/Lipnikov/Simonicini M3AS 2005
		// (1) orthonormalize columns of the matrix R
		Scalar norm = R[0][0]*R[0][0];
		for (int i = 1; i < sfs.size(); i++)
			norm += R[i][0]*R[i][0];
		norm = sqrt(norm);
		for (int i = 0; i < sfs.size(); i++)
			R[i][0] /= norm;
		Scalar weight = R[0][1]*R[0][0];
		for (int i = 1; i < sfs.size(); i++)
			weight += R[i][1]*R[i][0];
		for (int i = 0; i < sfs.size(); i++)
			R[i][1] -= weight*R[i][0];
		norm = R[0][1]*R[0][1];
		for (int i = 1; i < sfs.size(); i++)
			norm += R[i][1]*R[i][1];
		norm = sqrt(norm);
		for (int i = 0; i < sfs.size(); i++)
			R[i][1] /= norm;
		if (dim == 3)
		{
			Scalar weight1 = R[0][2]*R[0][0];
			Scalar weight2 = R[0][2]*R[0][1];
			for (int i = 1; i < sfs.size(); i++)
			{
				weight1 += R[i][2]*R[i][0];
				weight2 += R[i][2]*R[i][1];
			}
			for (int i = 0; i < sfs.size(); i++)
				R[i][1] -= weight1*R[i][0] + weight2*R[i][1];
			norm = R[0][2]*R[0][2];
			for (int i = 1; i < sfs.size(); i++)
				norm += R[i][2]*R[i][2];
			norm = sqrt(norm);
			for (int i = 0; i < sfs.size(); i++)
				R[i][2] /= norm;
		}
		//      std::cout << "~R =\dim" << R;

		// (2) Build the matrix ~D
		FieldMatrix<Scalar,2*dim,2*dim> D(0);
		for (int s = 0; s < sfs.size(); s++)
		{
			Dune::FieldVector<Scalar,2*dim> es(0);
			es[s] = 1;
			for (int k = 0; k < sfs.size(); k++)
			{
				D[k][s] = es[k];
				for (int i = 0; i < dim; i++)
				{
					D[k][s] -= R[s][i]*R[k][i];
				}
			}
		}

		Scalar traceK = K[0][0];
		for (int i = 1; i < dim; i++)
			traceK += K[i][i];
		D *= 2.0*traceK/volume;
		//      std::cout << "u~D =\dim" << D;

		// (3) Build the matrix W = Minv
		FieldMatrix<Scalar,2*dim,dim> NK(N);
		NK.rightmultiply (K);
		for (int i = 0; i < sfs.size(); i++)
		{
			for (int j = 0; j < sfs.size(); j++)
			{
				W[i][j] = NK[i][0]*N[j][0];
				for (int k = 1; k < dim; k++)
					W[i][j] += NK[i][k]*N[j][k];
			}
		}

		W /= volume;
		W += D;
		//       std::cout << "W = \dim" << W;
		//       std::cout << D[2][2] << ", " << D[2][0] << ", " << D[2][3] << ", " << D[2][1] << std::endl;
		//       std::cout << D[0][2] << ", " << D[0][0] << ", " << D[0][3] << ", " << D[0][1] << std::endl;
		//       std::cout << D[3][2] << ", " << D[3][0] << ", " << D[3][3] << ", " << D[3][1] << std::endl;
		//       std::cout << D[1][2] << ", " << D[1][0] << ", " << D[1][3] << ", " << D[1][1] << std::endl;


		// Now the notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
		// The matrix W developed so far corresponds to one element-associated
		// block of the matrix B^{-1} there.

		// Corresponding to the element under consideration,
		// calculate the part of the matrix C coupling velocities and element pressures.
		// This is just a row vector of size sfs.size().
		// scale with volume
		for (int i = 0; i < sfs.size(); i++)
			c[i] = faceVol[i];

		// Set up the element part of the matrix \Pi coupling velocities
		// and pressure-traces. This is a diagonal matrix with entries given by faceVol.
		for (int i = 0; i < sfs.size(); i++)
			Pi[i][i] = faceVol[i];

		// Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
		Dune::FieldVector<Scalar,2*dim> Wc(0);
		W.umv(c, Wc);
		dinv = 1.0/(c*Wc);

		// Calculate the element part of the matrix F = Pi W c^T which is a column vector.
		F = 0;
		Pi.umv(Wc, F);
		//      std::cout << "Pi = \dim" << Pi << "c = " << c << ", F = " << F << std::endl;

		// Calculate the source f
		int p = 0;
		qmean = 0;
		for (size_t g=0; g<Dune::QuadratureRules<Scalar,dim>::rule(gt,p).size(); ++g) // run through all quadrature points

		{
			const Dune::FieldVector<Scalar,dim>& local = Dune::QuadratureRules<Scalar,dim>::rule(gt,p)[g].position(); // pos of integration point
			Dune::FieldVector<Scalar,dim> global = element.geometry().global(local); // ip in global coordinates
			double weight = Dune::QuadratureRules<Scalar,dim>::rule(gt,p)[g].weight();// weight of quadrature point
			Scalar detjac = element.geometry().integrationElement(local); // determinant of jacobian
			Scalar factor = weight*detjac;
			Scalar q = problem_.sourcePress(global,element,local);
			qmean += q*factor;
		}

	}

private:
	void assembleV (const Element& element, int k=1)
	{
		// extract some important parameters
		Dune::GeometryType gt = element.geometry().type();
		const typename Dune::CRShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
		sfs=Dune::CRShapeFunctions<Scalar,Scalar,dim>::general(gt,1);
		setcurrentsize(sfs.size());

		// The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
		// The matrix W developed here corresponds to one element-associated
		// block of the matrix B^{-1} there.
		Dune::FieldVector<Scalar,2*dim> faceVol(0);
		Dune::FieldMatrix<Scalar,2*dim,2*dim> W(0);
		Dune::FieldVector<Scalar,2*dim> c(0);
		Dune::FieldMatrix<Scalar,2*dim,2*dim> Pi(0);
		Dune::FieldVector<Scalar,2*dim> F(0);
		Scalar dinv;
		Scalar qmean;
		this->assembleElementMatrices(element, faceVol, W, c, Pi, dinv, F, qmean);

		// Calculate the element part of the matrix Pi W Pi^T.
		Dune::FieldMatrix<Scalar,2*dim,2*dim> PiWPiT(W);
		PiWPiT.rightmultiply(Pi);
		PiWPiT.leftmultiply(Pi);

		// Calculate the element part of the matrix F D^{-1} F^T.
		Dune::FieldMatrix<Scalar,2*dim,2*dim> FDinvFT(0);
		for (int i = 0; i < sfs.size(); i++)
			for (int j = 0; j < sfs.size(); j++)
				FDinvFT[i][j] = dinv*F[i]*F[j];

		// Calculate the element part of the matrix S = Pi W Pi^T - F D^{-1} F^T.
		for (int i = 0; i < sfs.size(); i++)
			for (int j = 0; j < sfs.size(); j++)
				this->A[i][j] = PiWPiT[i][j] - FDinvFT[i][j];

		// Calculate the source term F D^{-1} f
		// NOT WORKING AT THE MOMENT
		Scalar factor = dinv*qmean;
		for (int i = 0; i < sfs.size(); i++)
			this->b[i] = F[i]*factor;

		//        std::cout << "faceVol = " << faceVol << std::endl << "W = " << std::endl << W << std::endl
		//              << "c = " << c << std::endl << "Pi = " << std::endl << Pi << std::endl
		//              << "dinv = " << dinv << std::endl << "F = " << F << std::endl;
		//        std::cout << "dinvF = " << dinvF << ", q = " << qmean
		//             << ", b = " << this->b[0] << ", " << this->b[1] << ", " << this->b[2] << ", " << this->b[3] << std::endl;
	}

	void assembleBC (const Element& element, int k=1)
	{
		// extract some important parameters
		Dune::GeometryType gt = element.geometry().type();
		const typename Dune::CRShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
		sfs=Dune::CRShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
		setcurrentsize(sfs.size());

		// determine quadrature order
		int p=0;

		// evaluate boundary conditions via intersection iterator
		typedef typename GridView::IntersectionIterator IntersectionIterator;

		//std::cout << "new element." << std::endl;
		IntersectionIterator endit = gridView_.template iend(element);
		for (IntersectionIterator it = gridView_.template ibegin(element); it!=endit; ++it)
		{
			//std::cout << "\tnew intersection iterator." << std::endl;

			// if we have a neighbor then we assume there is no boundary (forget interior boundaries)
			// in level assemble treat non-level neighbors as boundary
			if (it->neighbor())
			{
				if (levelBoundaryAsDirichlet_ && it->outside()->level()==element.level())
					continue;
				if (!levelBoundaryAsDirichlet_)
					continue;
			}

			//std::cout << "\t\tsurvived first if statements" << std::endl;

			// determine boundary condition type for this face, initialize with processor boundary
			typename BoundaryConditions::Flags bctypeface = BoundaryConditions::process;

			// handle face on exterior boundary, this assumes there are no interior boundaries
			//if (it.boundary())
			if (!it->neighbor())
			{
				//std::cout << "\t\t\tsurvived second if-statements." << std::endl;
				Dune::GeometryType gtface = it->geometryInInside().type();
				for (size_t g = 0; g < Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p).size(); ++g)
				{
					const Dune::FieldVector<Scalar,dim-1>& faceLocalNm1 = Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p)[g].position();
					FieldVector<Scalar,dim> local = it->geometryInInside().global(faceLocalNm1);
					FieldVector<Scalar,dim> global = it->geometry().global(faceLocalNm1);
					bctypeface = problem_.bctypePress(global,element,local); // eval bctype
					//std::cout << "\t\t\tlocal = " << local << ", global = " << global << ", bctypeface = " << bctypeface
					//        << ", size = " << Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p).size() << std::endl;


					if (bctypeface!=BoundaryConditions::neumann) break;

					Scalar J = problem_.neumannPress(global,element,local);
					double weightface = Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p)[g].weight();
					Scalar detjacface = it->geometry().integrationElement(faceLocalNm1);
					for (int i=0; i<sfs.size(); i++) // loop over test function number
						if (this->bctype[i][0]==BoundaryConditions::neumann)
						{
							this->b[i] -= J*sfs[i].evaluateFunction(0,local)*weightface*detjacface;
						}
				}
				if (bctypeface==BoundaryConditions::neumann) continue; // was a neumann face, go to next face
			}

			// If we are here, then it is
			// (i)   an exterior boundary face with Dirichlet condition, or
			// (ii)  a processor boundary (i.element. neither boundary() nor neighbor() was true), or
			// (iii) a level boundary in case of level-wise assemble
			// How processor boundaries are handled depends on the processor boundary mode
			if (bctypeface==BoundaryConditions::process && procBoundaryAsDirichlet_==false
					&& levelBoundaryAsDirichlet_==false)
				continue; // then it acts like homogeneous Neumann

			// now handle exterior or interior Dirichlet boundary
			for (int i=0; i<sfs.size(); i++) // loop over test function number

			{
				if (sfs[i].codim()==0) continue; // skip interior dof
				if (sfs[i].codim()==1) // handle face dofs

				{
					if (sfs[i].entity()==it->indexInInside())
					{
						if (this->bctype[i][0]<bctypeface)
						{
							this->bctype[i][0] = bctypeface;
							if (bctypeface==BoundaryConditions::process)
								this->b[i] = 0;
							if (bctypeface==BoundaryConditions::dirichlet)
							{
								Dune::FieldVector<Scalar,dim> global = element.geometry().global(sfs[i].position());
								//std::cout << "i = " << i << ", loop 1, global = " << global << std::endl;
								this->b[i] = problem_.dirichletPress(global,element,sfs[i].position());
							}
						}
					}
					continue;
				}
				// handle subentities of this face
				for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it->indexInInside(),1,sfs[i].codim()); j++)
					if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it->indexInInside(),1,j,sfs[i].codim()))
					{
						if (this->bctype[i][0]<bctypeface)
						{
							this->bctype[i][0] = bctypeface;
							if (bctypeface==BoundaryConditions::process)
								this->b[i] = 0;
							if (bctypeface==BoundaryConditions::dirichlet)
							{
								Dune::FieldVector<Scalar,dim> global = element.geometry().global(sfs[i].position());
								//std::cout << "loop 2, global = " << global << std::endl;
								this->b[i] = problem_.dirichletPress(global,element,sfs[i].position());
							}
						}
					}
			}
		}
	}

	// parameters given in constructor

private:
	Problem& problem_;
	bool levelBoundaryAsDirichlet_;
	bool procBoundaryAsDirichlet_;
	const GridView& gridView_;
};

/** @} */
}
#endif
