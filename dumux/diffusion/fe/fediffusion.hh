// $Id$ 

#ifndef DUNE_FEDIFFUSION_HH
#define DUNE_FEDIFFUSION_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/diffusion/diffusion_deprecated.hh"
#include "dune/disc/operators/p1operator.hh"
#include "dumux/diffusion/fe/p1groundwater.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
//#include "mgpreconditioner.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 */

namespace Dune
{
  //! \ingroup diffusion
  //! Base class for defining an instance of a numerical diffusion model.
  /*! An interface for defining a numerical diffusion model for the
   *  solution of equations of the form
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
   * \f$p = g\f$ on \f$\Gamma_1\f$, and
   * \f$\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$
   * on \f$\Gamma_2\f$. Here,
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
   * and \f$\lambda\f$ the total mobility, possibly depending on the
   * saturation, \f$q\f$ the source term.
	Template parameters are:

	- Grid      a DUNE grid type
	- RT        type used for return values
   */
  template<class G, class RT, class VC, class LocalStiffnessType = GroundwaterEquationLocalStiffness<G,RT> >
  class FEDiffusion
  : public Diffusion< G, RT, VC >
  {
	  template<int dim>
	  struct ElementLayout
	  {
		  bool contains (GeometryType gt)
	      {
			  return gt.dim() == dim;
	      }
	  };

	  typedef LevelP1Function<G,RT,1> PressP1Type;
	  typedef LevelP1OperatorAssembler<G,RT,1> LevelOperatorAssembler;

  public:
        typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > VelType;
	typedef BlockVector< FieldVector<RT,1> > RepresentationType;

	void assemble(const RT t=0)
	{
		LocalStiffnessType lstiff(this->problem, false, this->grid, this->level());
		A.assemble(lstiff, pressP1, f);
		return;
	}

	void solve()
	{
		typedef typename LevelP1Function<G,RT>::RepresentationType VectorType;
		typedef typename LevelP1OperatorAssembler<G,RT,1>::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;

		Operator op(*A);  // make operator out of matrix
		double red=1E-10;
		SeqILU0<MatrixType,VectorType,VectorType> ilu0(*A,1.0);// a precondtioner
		CGSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
		InverseOperatorResult r;
		solver.apply(*pressP1, *f, r);
		this->press = *pressP1;

		return;
	}

	void pressure(const RT t=0)
	{
		assemble(t);
		solve();
		return;
	}

	void totalVelocity(VelType& velocity, const RT t=0) const
	{
		DUNE_THROW(NotImplemented,"FEDiffusion :: totalVelocity");
		return;
	}

	void vtkout (const char* name, int k) const
	{
		VTKWriter<typename G::LevelGridView>
			vtkwriter(this->grid.levelView(this->level()));
		char fname[128];
		sprintf(fname,"%s-%05d",name,k);
		vtkwriter.addVertexData(this->diffproblem.variables.pressure,"total pressure p~");
		vtkwriter.write(fname, VTKOptions::ascii);
	}

	FEDiffusion(G& g, DiffusionProblem<G, RT, VC>& prob,
		    TransportProblem<G, RT, VC>& satprob = *(new typename Dune::SimpleProblem<G, RT, VC>),
		    int lev = -1)
	  : Diffusion<G, RT, VC>(g, prob, lev == -1 ? g.maxLevel() : lev),
	  pressP1(g, this->level()), f(g, this->level()), A(g, this->level())
	{
		this->press.resize(g.size(this->level(), G::dimension));
		*pressP1 = 0;
	}

	PressP1Type pressP1;
	PressP1Type f;
	LevelOperatorAssembler A;
  };

}
#endif
