#ifndef DUNE_2P2CBOXJACOBIAN_HH
#define DUNE_2P2CBOXJACOBIAN_HH

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
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/functions/p1function.hh>

#include"localjacobian_2p2c.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */

namespace Dune {
/** @addtogroup DISC_Disc
 *
 * @{
 */
/**
 * @brief compute local jacobian matrix for conforming finite elements for diffusion equation
 *
 */

//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the 
 diffusion equation

 div j = q; j = -K grad u; in Omega

 u = g on Gamma1; j*n = J on Gamma2.

 Uses conforming finite elements with the Lagrange shape functions.
 It should work for all dimensions and element types.
 All the numbering is with respect to the reference element and the
 Lagrange shape functions

 Template parameters are:

 - Grid  a DUNE grid type
 - RT    type used for return values 
 */
template<class Imp, class G, class RT, int m,
		class BoxFunction = LeafP1Function<G, RT, m> > class BoxJacobian2p2c :
	public LocalJacobian2p2c<Imp,G,RT,m> {
	// mapper: one data element per vertex
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == 0;
		}
	};

	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename Entity::Geometry Geometry;
	typedef typename LocalJacobian2p2c<Imp,G,RT,m>::VBlockType VBlockType;
	typedef typename LocalJacobian2p2c<Imp,G,RT,m>::MBlockType MBlockType;
	typedef FVElementGeometry<G> FVElementGeometry;
	typedef MultipleCodimMultipleGeomTypeMapper<G, typename G::Traits::LeafIndexSet, P1Layout>
			VertexMapper;

public:
	// define the number of components of your system, this is used outside
	// to allocate the correct size of (dense) blocks with a FieldMatrix
	enum {n=G::dimension};
	enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};

	//! Constructor
	BoxJacobian2p2c(bool levelBoundaryAsDirichlet_, const G& grid,
			BoxFunction& sol, bool procBoundaryAsDirichlet_=true) :
		vertexMapper(grid, grid.leafIndexSet()),
				levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
				procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
				currentSolution(sol), oldSolution(grid), dt(1) {
	}

	//**********************************************************
	//*																			*
	//*	Computation of the local defect								*
	//*																			*
	//**********************************************************

	template<class TypeTag> void localDefect(const Entity& e,
			const VBlockType* sol) {
		for (int i=0; i < this->fvGeom.nNodes; i++)
			this->def[i] = 0;

		computeElementData(e);
		updateVariableData(e, sol);

		for (int i=0; i < this->fvGeom.nNodes; i++) // begin loop over vertices / sub control volumes
		{
			// implicit Euler
			VBlockType massContrib = computeM(e, sol, i);
			massContrib -= computeM(e, uold, i);
			massContrib *= this->fvGeom.subContVol[i].volume/dt;
			this->def[i] += massContrib;
			std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
			std::cout.setf(std::ios_base::uppercase);
			std::cout.precision(3);
			//std::cout << "i = " << i << ", massContrib = " << massContrib << std::endl;

			// get source term 
			VBlockType q = computeQ(e, sol, i);
			q *= this->fvGeom.subContVol[i].volume;
			this->def[i] -= q;
		} // end loop over vertices / sub control volumes

		for (int k = 0; k < this->fvGeom.nEdges; k++) // begin loop over edges / sub control volume faces
		{
			int i = this->fvGeom.subContVolFace[k].i;
			int j = this->fvGeom.subContVolFace[k].j;

			VBlockType flux = computeA(e, sol, k);

			// add to defect 
			this->def[i] -= flux;
			this->def[j] += flux;
			//std::cout << "i = " << i << ", j = " << j << ", flux = " << flux << std::endl;
		} // end loop over edges / sub control volume faces

		// assemble boundary conditions 
		assembleBC<TypeTag> (e);

		// add to defect 
		for (int i=0; i < this->fvGeom.nNodes; i++) {
			this->def[i] += this->b[i];
//			 		  std::cout << ", b[" << i << "] = " << this->b[i] << std::endl;
		}
		// 		std::cout << std::endl;

		return;
	}

	void setLocalSolution(const Entity& e) {
		Dune::GeometryType gt = e.geometry().type();
		const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type
				&sfs=Dune::LagrangeShapeFunctions<DT, RT, n>::general(gt, 1);
		int size = sfs.size();
		this->setcurrentsize(size);

		for (int i = 0; i < size; i++)
			for (int comp = 0; comp < m; comp++) {
				this->u[i][comp]= currentSolution.evallocal(comp, e,
						sfs[i].position());
				uold[i][comp]
						= oldSolution.evallocal(comp, e, sfs[i].position());
			}

		return;
	}

	void setDt(double d) {
		dt = d;

		return;
	}

	double getDt() {
		return dt;
	}

	void setOldSolution(BoxFunction& uOld) {
		*oldSolution = *uOld;
	}

	VBlockType computeM(const Entity& e, const VBlockType* sol, int node) {
		return this->getImp().computeM(e, sol, node);
	}

	VBlockType computeQ(const Entity& e, const VBlockType* sol, int node) {
		return this->getImp().computeQ(e, sol, node);
	}

	VBlockType computeA(const Entity& e, const VBlockType* sol, int face) {
		return this->getImp().computeA(e, sol, face);
	}

	void computeElementData(const Entity& e) {
		return this->getImp().computeElementData(e);
	}

	// analog to EvalStaticData in MUFTE
	virtual void updateStaticData(const Entity& e, const VBlockType* sol) {
		return this->getImp().updateStaticData(e, sol);
	}

	// analog to EvalPrimaryData in MUFTE, uses members of varNData
	virtual void updateVariableData(const Entity& e, const VBlockType* sol) {
		return this->getImp().updateVariableData(e, sol);
	}

	template<class TypeTag> void assembleBC(const Entity& e) {
		Dune::GeometryType gt = e.geometry().type();
		const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type
				&sfs=Dune::LagrangeShapeFunctions<DT, RT, n>::general(gt, 1);
		setcurrentsize(sfs.size());
		this->fvGeom.update(e);

		const typename ReferenceElementContainer<DT,n>::value_type
				&referenceElement = ReferenceElements<DT, n>::general(gt);

		for (int i = 0; i < sfs.size(); i++) {
			this->bctype[i].assign(BoundaryConditions::neumann);
			this->b[i] = 0;
			this->dirichletIndex[i] = 0;
		}

		// evaluate boundary conditions via intersection iterator
		typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator
				IntersectionIterator;

		IntersectionIterator
				endit = IntersectionIteratorGetter<G, TypeTag>::end(e);
		for (IntersectionIterator
				it = IntersectionIteratorGetter<G, TypeTag>::begin(e); it
				!=endit; ++it) {
			// if we have a neighbor then we assume there is no boundary (forget interior boundaries)
			// in level assemble treat non-level neighbors as boundary
			if (it->neighbor()) {
				if (levelBoundaryAsDirichlet && it->outside()->level()==e.level())
					continue;
				if (!levelBoundaryAsDirichlet)
					continue;
			}

			// determine boundary condition type for this face, initialize with processor boundary
			FieldVector<typename BoundaryConditions::Flags, m> bctypeface(BoundaryConditions::process);
			FieldVector<int,m> dirichletIdx;

			// handle face on exterior boundary, this assumes there are no interior boundaries
			if (it->boundary()) {
				int faceIdx = it->numberInSelf();
				// 				std::cout << "faceIdx = " << faceIdx << ", beginning: " << std::endl;
				// 				for (int i = 0; i < 4; i++) 
				// 				  std::cout << "bctype[" << i << "] = " << this->bctype[i] << std::endl; 

				int nNodesOfFace = referenceElement.size(faceIdx, 1, n);
				for (int nodeInFace = 0; nodeInFace < nNodesOfFace; nodeInFace++) {
					int nodeInElement = referenceElement.subEntity(faceIdx, 1,
							nodeInFace, n);
					for (int equationNumber = 0; equationNumber < m; equationNumber++) {
						if (this->bctype[nodeInElement][equationNumber]
								== BoundaryConditions::neumann) {
							int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,
									nodeInFace);
							FieldVector<DT,n>
									local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
							FieldVector<DT,n>
									global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
							bctypeface = this->getImp().problem.bctype(global, e, it, local); // eval bctype
							this->getImp().problem.dirichletIndex(global, e, it, local, dirichletIdx); // eval bctype
//							 						std::cout << "faceIdx = " << faceIdx << ", nodeInElement = " << nodeInElement 
//							 							  << ", bfIdx = " << bfIdx << ", local = " << local << ", global = " << global 
//							 							  << ", bctypeface = " << bctypeface << std::endl; 
							if (bctypeface[equationNumber]
									!=BoundaryConditions::neumann)
								break;
							VBlockType J = this->getImp().problem.J(global, e, it, local);
							J[equationNumber]
									*= this->fvGeom.boundaryFace[bfIdx].area;
							this->b[nodeInElement][equationNumber]
									+= J[equationNumber];
						}
					}
				}

				if (bctypeface[0]==BoundaryConditions::neumann && bctypeface[1]
						==BoundaryConditions::neumann)
					continue; // was a neumann face, go to next face
			}

			// If we are here, then it is 
			// (i)   an exterior boundary face with Dirichlet condition, or
			// (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
			// (iii) a level boundary in case of level-wise assemble
			// How processor boundaries are handled depends on the processor boundary mode
			if (bctypeface[0]==BoundaryConditions::process
					&& procBoundaryAsDirichlet==false
					&& levelBoundaryAsDirichlet==false)
				continue; // then it acts like homogeneous Neumann
			if (bctypeface[1]==BoundaryConditions::process
					&& procBoundaryAsDirichlet==false
					&& levelBoundaryAsDirichlet==false)
				continue; // then it acts like homogeneous Neumann

			for (int equationNumber=0; equationNumber<m; equationNumber++) {
				for (int i=0; i<sfs.size(); i++) // loop over test function number
				{
					//std::cout<<"i = "<<i<<std::endl;
					if (sfs[i].codim()==0)
						continue; // skip interior dof
					if (sfs[i].codim()==1) // handle face dofs
					{
						if (sfs[i].entity()==it->numberInSelf()) {
							if (this->bctype[i][equationNumber]
									<bctypeface[equationNumber]) {
								this->bctype[i][equationNumber]
										= bctypeface[equationNumber];
								this->dirichletIndex[i][equationNumber]
										= dirichletIdx[equationNumber];
								if (bctypeface[equationNumber]
										==BoundaryConditions::process)
									this->b[i][equationNumber] = 0;
								if (bctypeface[equationNumber]
										==BoundaryConditions::dirichlet) {
									this->b[i][equationNumber] = 0;
								}
							}
						}
						continue;
					}
					// handle subentities of this face
					for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it->numberInSelf(), 1, sfs[i].codim()); j++)
						if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it->numberInSelf(), 1, j,
								sfs[i].codim())) {
							if (this->bctype[i][equationNumber]
									<bctypeface[equationNumber]) {
								this->bctype[i][equationNumber]
										= bctypeface[equationNumber];
								this->dirichletIndex[i][equationNumber]
										= dirichletIdx[equationNumber];
								if (bctypeface[equationNumber]
										==BoundaryConditions::process)
									this->b[i][equationNumber] = 0;
								if (bctypeface[equationNumber]
										==BoundaryConditions::dirichlet) {
									this->b[i][equationNumber] = 0;
								}
							}
						}
				}
			}
		}

	}

	// parameters given in constructor
	VertexMapper vertexMapper;
	bool levelBoundaryAsDirichlet;
	bool procBoundaryAsDirichlet;
	const BoxFunction& currentSolution;
	BoxFunction oldSolution;

public:
	double dt;
	VBlockType uold[SIZE];
};
}
#endif
