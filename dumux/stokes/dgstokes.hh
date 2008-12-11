// $Id$
#ifndef DUNE_DGSTOKES_HH
#define DUNE_DGSTOKES_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/io.hh>

#include<dune/grid/common/quadraturerules.hh>

#include<dune/disc/shapefunctions/dgspace/monomialshapefunctions.hh>

#include"stokesparameters.hh"
#include"stokesproblem.hh"
#include"boundaryconditions.hh"
#include"testfunctions.hh"

#include <dune/istl/solvers.hh>
#include "dumux/pardiso/pardiso.hh"

namespace Dune
{
template<class G,int v_order, int p_order>
  class DGFiniteElementMethod
  {
	//dimension of grid
	enum {dim=G::dimension};
	enum { dimw=G::dimensionworld };
  public:
	//Grid
	typedef G Grid;
	//coordinate type
    typedef typename Grid::ctype ctype;
	typedef Dune::FieldVector< double , dim> Gradient;
	typedef Dune::FieldMatrix< double , dim, dim> InverseJacobianMatrix;
	// "order" is order of velocity shapefn
	// order of pressure shapefn  is (order-1) i.e, assumed one order less than that of velocity
	enum{v_ordr = v_order};
	enum{p_ordr=p_order};

	//local vector and matrix blocks
	static const int BlockSize =dim*Dune::MonomialShapeFunctionSetSize<dim,v_order>::maxSize+Dune::MonomialShapeFunctionSetSize<dim,p_order>::maxSize;
	static const int VBlockSize =((dim*Dune::MonomialShapeFunctionSetSize<dim,v_order>::maxSize));
	static const int PBlockSize =((Dune::MonomialShapeFunctionSetSize<dim,p_order>::maxSize));

	typedef Dune::FieldVector<double,BlockSize> LocalVectorBlock;
	typedef Dune::FieldMatrix<double,BlockSize,BlockSize> LocalMatrixBlock;
	//shapefn
	typedef Dune::MonomialShapeFunctionSet<ctype,double,dim> ShapeFunctionSet;
	inline const ShapeFunctionSet & getVelocityShapeFunctionSet(Dune::GeometryType gt) const;
	inline const ShapeFunctionSet & getPressureShapeFunctionSet(Dune::GeometryType gt) const;

	typedef typename Grid::template Codim<0>::LeafIterator ElementLeafIterator;
	typedef typename Grid::template Codim<0>::LevelIterator ElementLevelIterator;
	typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
	typedef typename Grid::template Codim<0>::Entity Entity;
	typedef typename Grid::template Codim<dim>::LevelIterator VertexIterator;
	typedef typename Grid::template Codim<dim>::Entity Vertex;
	typedef typename Grid::template Codim<0>::LevelIntersectionIterator IntersectionLevelIterator;
	typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
	typedef typename Grid::template Codim<1>::EntityPointer InterSectionPointer;




    DGFiniteElementMethod (Grid &g, StokesProblem<Grid, ctype>& prob, const DGStokesParameters& par)
    : grid(g), problem_(prob), parameter(par)
    {}

	//local assembly
	void assembleVolumeTerm(Entity& ep, LocalMatrixBlock& Aee,LocalVectorBlock& Be) const;
	void assembleFaceTerm(Entity& ep,IntersectionIterator& isp, LocalMatrixBlock& Aee,
						  LocalMatrixBlock& Aef,LocalMatrixBlock& Afe, LocalVectorBlock& Be) const;

	void assembleDirichletBoundaryTerm(Entity& ep, IntersectionIterator& isp, LocalMatrixBlock& Aee,LocalVectorBlock& Be)const ;
	void assembleNeumannBoundaryTerm(Entity& ep, IntersectionIterator& isp, LocalMatrixBlock& Aee,LocalVectorBlock& Be)const ;
	void assembleInterfaceTerm(Entity& ep, IntersectionIterator& isp, LocalMatrixBlock& Aee,LocalVectorBlock& Be)const ;
	// stokes system has dim+1 variables (dim velocity comps and 1 pressure)
	double evaluateSolution(int variable,const Entity& element,const Dune::FieldVector<ctype,dim>& local,
							const LocalVectorBlock& xe) const;
	Gradient evaluateGradient(int variable,const Entity& element,const Dune::FieldVector<ctype,dim>& local,
							const LocalVectorBlock& xe) const;
	double evaluateL2error(int variable, const Entity& element, const LocalVectorBlock& xe)const;

	double evaluateH1error(int variable, const Entity& element, const LocalVectorBlock& xe)const;

    StokesProblem<Grid, ctype>& problem()
    {
      return problem_;
    }


  private:
	Grid& grid;
    StokesProblem<Grid, ctype>& problem_;
	Dune::MonomialShapeFunctionSetContainer<ctype,double,dim,v_order> vspace;
	Dune::MonomialShapeFunctionSetContainer<ctype,double,dim,(p_order)> pspace;
	DGStokesParameters parameter;
  };



template<class G,int v_order,int p_order>
  class DGStokes
  {

  public:
	//dimension of grid
	enum {dimension=G::dimension};
	enum { dimensionworld=G::dimensionworld };
	enum{v_ordr = v_order};
	enum{p_ordr=p_order};
	typedef G Grid;
  private:
	enum {dim=G::dimension};
	enum { dimw=G::dimensionworld };
    typedef typename Grid::ctype ctype;
	// Iterators
	typedef typename Grid::template Codim<0>::LeafIterator ElementLeafIterator;
	typedef typename Grid::template Codim<0>::LevelIterator ElementLevelIterator;
	typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
	typedef typename Grid::template Codim<0>::Entity Entity;
	typedef typename Grid::template Codim<dim>::Entity Vertex;
	typedef typename Grid::template Codim<dim>::LevelIterator VertexIterator;
	typedef typename Grid::template Codim<0>::LevelIntersectionIterator IntersectionLevelIterator;
	typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
	typedef typename Grid::template Codim<1>::EntityPointer InterSectionPointer;
	static const int BlockSize =((dim*Dune::MonomialShapeFunctionSetSize<dim,v_order>::maxSize)+(Dune::MonomialShapeFunctionSetSize<dim,p_order>::maxSize));

	static const int VBlockSize =((dim*Dune::MonomialShapeFunctionSetSize<dim,v_order>::maxSize));
	static const int PBlockSize =((Dune::MonomialShapeFunctionSetSize<dim,p_order>::maxSize));

	typedef typename DGFiniteElementMethod<G,v_order,p_order>::Gradient Gradient;
	typedef typename DGFiniteElementMethod<G,v_order,p_order>::LocalVectorBlock LocalVectorBlock;
	typedef typename DGFiniteElementMethod<G,v_order,p_order>::LocalMatrixBlock LocalMatrixBlock;
	typedef Dune::BlockVector<LocalVectorBlock> LVector;
	typedef Dune::BCRSMatrix<LocalMatrixBlock> LMatrix;

  public:
	typedef LMatrix MatrixType;
	typedef G GridType;

	  MatrixType& matrix()
	  {
		  return A;
	  }

	  LVector& rhs()
	  {
		  return b;
	  }

	  LVector& sol()
	  {
		  return solution;
	  }

	  void initial()
	  {}

    DGStokes(Grid &g, StokesProblem<Grid, ctype>& prob, DGStokesParameters& par)
      : grid(g),level(g.maxLevel()), dgfem(g,prob,par)
	{
	  if (par.sigma==1 & par.epsilon==1)
	std::cout<<"You are using NIPG scheme"<<std::endl;
  else if(par.sigma==1 & par.epsilon==-1)
	std::cout<<"You are using SIPG scheme"<<std::endl;
  else if(par.sigma==0 & par.epsilon==1)
	std::cout<<"You are using OBB scheme"<<std::endl;
  else
	std::cout<<"check DG parameters epsilon and sigma"<<std::endl;
	};

	// global assembly and solving
	void assembleStokesSystem() ;
	void assemble()
	{
		assembleStokesSystem();
	}
	void solveStokesSystem();
	double evaluateSolution(const EntityPointer & e,
                            const Dune::FieldVector<ctype, dim> & local) const;

	//l2error computation
	// stokes system has dim+1 variables (dim velocity comps and 1 pressure)
	double l2errorStokesSystem(int variable) const;
	double h1errorStokesSystem(int variable) const;
	void convertToCellData(int variable, BlockVector<FieldVector<double, 1> >& cellData);
	void vtkout (const char* name, int k);

  private:
	typedef typename DGFiniteElementMethod<G,v_order,p_order>::ShapeFunctionSet ShapeFunctionSet;
	inline const ShapeFunctionSet & getVelocityShapeFunctionSet(const EntityPointer &ep) const;
	inline const ShapeFunctionSet & getPressureShapeFunctionSet(const EntityPointer &ep) const;
  public:
	Grid & grid;
	int level;
	DGFiniteElementMethod<G,v_order,p_order> dgfem;
	LMatrix A;
    LVector b;
	LVector solution;
  };



#include "dgstokes.cc"

} // end of namespace

#endif
