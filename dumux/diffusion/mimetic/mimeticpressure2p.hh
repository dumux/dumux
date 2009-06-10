// $Id$

#ifndef DUNE_MIMETICPRESSURE2P_HH
#define DUNE_MIMETICPRESSURE2P_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/diffusion/diffusion.hh"
#include "dumux/operators/mimeticoperator.hh"
#include "dumux/diffusion/mimetic/mimeticgroundwater.hh"

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
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 Template parameters are:

 - Grid      a DUNE grid type
 - RT        type used for return values
 */
template<
        class GridView,
        class Scalar,
        class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = Dune::MimeticGroundwaterEquationLocalStiffness<GridView,Scalar,VC, Problem> >
class MimeticPressure2P: public Diffusion<GridView, Scalar, VC, Problem>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        Sw = 0, Sn = 1
    };
typedef    typename GridView::Grid Grid;
    typedef Dune::LevelCRFunction<Grid,Scalar,1> TraceType;
    typedef Dune::LevelP0Function<Grid,Scalar,2*GridView::dimension> NormalVelType;
    typedef Dune::MimeticOperatorAssembler<Grid,Scalar,1> LevelOperatorAssembler;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef BlockVector< Dune::FieldVector<Scalar,1> > RepresentationType;

    void assemble(const Scalar t=0)
    {
        LocalStiffnessType lstiff(this->diffProblem, false, this->gridView);
        A.assemble(lstiff, pressTrace, f);
        return;
    }

    void solve();

    void postprocess()
    {
        LocalStiffnessType lstiff(this->diffProblem, false, this->gridView);
        A.calculatePressure(lstiff, pressTrace, normalVelocity, this->diffProblem.variables().pressure());
        //printvector(std::cout, this->variables.pressure, "element pressures", "row", 200, 1, 5);
        //printvector(std::cout, *normalVelocity, "normal velocities", "row", 200, 1, 5);
        return;
    }

    void pressure(bool first = true, const Scalar t=0)
    {
        if (first)
        {
            initializeMaterialLaws();
        }
        assemble(t);
        solve();
        postprocess();
        return;
    }

    void calculateVelocity(const Scalar t) const;

    void calculateVelocity(const Scalar t, double lev) const
    {
        DUNE_THROW(Dune::NotImplemented, "upscaled velocities only implemented in FVDiffusion");
    }

    //constitutive functions are initialized and stored in the variables object
    void initializeMaterialLaws();

    void vtkout (const char* name, int k) const
    {
        this->diffProblem().variables().vtkout(name, k);
    }

    MimeticPressure2P(GridView& gridView, Problem& prob, std::string satType = "Sw", int level = -1, bool calcPressure = true)
    : Diffusion<GridView, Scalar, VC, Problem>(gridView, prob),
    saturationType((satType == "Sw") ? 0 : ((satType == "Sn") ? 1 : 999)),
    level_((level >= 0) ? level : gridView.grid().maxLevel()), pressTrace(gridView.grid(), level_), normalVelocity(gridView.grid(), level_), f(gridView.grid(), level_), A(gridView.grid(), level_)
    {
        *pressTrace = 0;
        *f = 0;
        if (saturationType == 999)
        {
            DUNE_THROW(NotImplemented, "Saturation type not supported!");
        }
    }

public:
    const int saturationType;
private:
    int level_;
public:
    TraceType pressTrace; //!< vector of pressure traces
    NormalVelType normalVelocity;
    TraceType f;
    LevelOperatorAssembler A;

};
template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType> void MimeticPressure2P<GridView, Scalar, VC, Problem,LocalStiffnessType>::solve()
{
    typedef typename LevelCRFunction<Grid,Scalar>::RepresentationType VectorType;
    typedef typename LevelCROperatorAssembler<Grid,Scalar,1>::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;

    //printmatrix(std::cout, *A, "global stiffness matrix", "row", 11, 3);
    //printvector(std::cout, *f, "right hand side", "row", 200, 1, 5);
    Operator op(*A); // make operator out of matrix
    double red=1E-12;
    SeqILU0<MatrixType,VectorType,VectorType> ilu0(*A,1.0);// a precondtioner
    //SeqJac<MatrixType,VectorType,VectorType> ilu0(*A,1,0.9);// a precondtioner
    //SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*A);// a precondtioner
    BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1); // an inverse operator
    //CGSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
    InverseOperatorResult r;
    solver.apply(*pressTrace, *f, r);
    //printvector(std::cout, *pressTrace, "solution", "row", 200, 1, 5);
    return;
}

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType> void MimeticPressure2P<GridView, Scalar, VC, Problem,LocalStiffnessType>::calculateVelocity(const Scalar t=0) const
{
    // ASSUMES axiparallel grids in 2D
    for (int i = 0; i < this->gridView.size(0); i++)
    {
        this->diffProblem.variables().velocity()[i][0][0] = -(*normalVelocity)[i][0];
        this->diffProblem.variables().velocity()[i][0][1] = 0;
        this->diffProblem.variables().velocity()[i][1][0] = (*normalVelocity)[i][1];
        this->diffProblem.variables().velocity()[i][1][1] = 0;
        this->diffProblem.variables().velocity()[i][2][0] = 0;
        this->diffProblem.variables().velocity()[i][2][1] = -(*normalVelocity)[i][2];
        this->diffProblem.variables().velocity()[i][3][0] = 0;
        this->diffProblem.variables().velocity()[i][3][1] = (*normalVelocity)[i][3];
    }
    return;
}

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType> void MimeticPressure2P<GridView, Scalar, VC, Problem,LocalStiffnessType>::initializeMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        int globalIdx = this->diffProblem.variables().indexDiffusion(*eIt);

        Scalar sat = this->diffProblem.variables().saturation()[globalIdx];

        std::vector<Scalar> mobilities(2);

        if (saturationType == Sw)
        {
            mobilities = this->diffProblem.materialLaw().mob(sat, globalPos, *eIt, localPos);
            this->diffProblem.variables().capillaryPressure()[globalIdx]= this->diffProblem.materialLaw().pC(sat, globalPos, *eIt, localPos);
        }
        else if (saturationType == Sn)
        {
            mobilities = this->diffProblem.materialLaw().mob(1-sat, globalPos, *eIt, localPos);
            this->diffProblem.variables().capillaryPressure()[globalIdx]= this->diffProblem.materialLaw().pC(1-sat, globalPos, *eIt, localPos);
        }
        else
        {
            DUNE_THROW(RangeError, "materialLaws not initialized!");
        }

        // initialize mobilities
        this->diffProblem.variables().mobilityWetting()[globalIdx]= mobilities[0];
        this->diffProblem.variables().mobilityNonWetting()[globalIdx]= mobilities[1];
        this->diffProblem.variables().fracFlowFuncWetting()[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
        this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
    }
    return;
}
}
#endif
