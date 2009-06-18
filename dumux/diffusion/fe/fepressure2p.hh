// $Id:$
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
#ifndef DUNE_FEPRESSURE2P_HH
#define DUNE_FEPRESSURE2P_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/disc/operators/p1operator.hh>

#include "dumux/diffusion/diffusion.hh"
#include "dumux/diffusion/fe/p1groundwater.hh"
#include "dumux/pardiso/pardiso.hh"

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
 - Scalar        type used for return values
 */
template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
class FEPressure2PBase: public Diffusion<GridView, Scalar, VC, Problem>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains(GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };
    typedef P1Function<GridView,Scalar,Communication,1> PressP1Type;
    typedef typename GridView::Grid Grid;
    typedef P1OperatorAssembler<Grid,Scalar,GridView,Communication,1> OperatorAssembler;

public:
    typedef BlockVector< FieldVector<FieldVector<Scalar, GridView::dimension>, 2*GridView::dimension> > VelType;
    typedef BlockVector< FieldVector<Scalar,1> > RepresentationType;

    void assemble(const Scalar t=0)
    {
        LocalStiffnessType lstiff(this->diffProblem, false, this->gridView);
        A.assemble(lstiff, pressP1, f);
        return;
    }

    void solve()
    {
    	std::cout << "FEPressure2P: solve for pressure" << std::endl;

        typedef typename PressP1Type::RepresentationType Vector;
        typedef typename OperatorAssembler::RepresentationType Matrix;
        typedef MatrixAdapter<Matrix,Vector,Vector> Operator;

        Operator op(*A);
        double reduction = 1E-12;
        int maxIt = 10000;
        int verboseLevel = 1;
        InverseOperatorResult result;

        if (preconditionerName_ == "SeqILU0")
        {
        	SeqILU0<Matrix,Vector,Vector> preconditioner(*A, 1.0);
        	if (solverName_ == "CG")
        	{
        		CGSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
        		solver.apply(*pressP1, *f, result);
        	}
        	else if (solverName_ == "BiCGSTAB")
        	{
        		BiCGSTABSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
        		solver.apply(*pressP1, *f, result);
        	}
        	else
        		DUNE_THROW(NotImplemented, "FEPressure2P :: solve : combination "
        				<< preconditionerName_<< " and "<< solverName_ << ".");
        }
        else if (preconditionerName_ == "SeqPardiso")
        {
        	SeqPardiso<Matrix,Vector,Vector> preconditioner(*A);
        	if (solverName_ == "Loop")
        	{
        		LoopSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
        		solver.apply(*pressP1, *f, result);
        	}
        	else
        		DUNE_THROW(NotImplemented, "FEPressure2P :: solve : combination "
        				<< preconditionerName_<< " and "<< solverName_ << ".");
        }
        else
        	DUNE_THROW(NotImplemented, "FEPressure2P :: solve : preconditioner "
        			<< preconditionerName_ << ".");

        this->diffProblem.variables().pressure() = *pressP1;

        return;
    }

    void pressure(const Scalar t=0)
    {
        assemble(t);
        solve();
        return;
    }

    void calculateVelocity(const Scalar t) const;

    void vtkout (const char* name, int k) const
    {
        this->diffProblem.variables().vtkout(name, k);
    }

    FEPressure2PBase(GridView& gridView, Problem& problem, Communication& comm, int level,
				std::string solver, std::string preconditioner)
    : Diffusion<GridView, Scalar, VC, Problem>(gridView, problem),
    level_((level >= 0) ? level : gridView.grid().maxLevel()),
    pressP1(gridView, comm), f(gridView, comm), A(gridView.grid(), gridView, comm),
    solverName_(solver), preconditionerName_(preconditioner)
    {
        *pressP1 = 0;
    }


private:
    int level_;
public:
    PressP1Type pressP1;
    PressP1Type f;
    OperatorAssembler A;
    std::string solverName_;
    std::string preconditionerName_;
};

template<class GridView, class Scalar, class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = GroundwaterEquationLocalStiffness<GridView, Scalar, Problem> >
class FEPressure2P: public FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid>  >
{
public:
    FEPressure2P(GridView& gridView, Problem& problem, int level = -1, std::string solver = "CG", std::string preconditioner = "SeqILU0")
    : FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >(
            gridView, problem,
            *(new LevelCommunicate<typename GridView::Grid>(gridView.grid(), level)),
            level, solver, preconditioner)
    {}

    FEPressure2P(GridView& gridView, Problem& problem, std::string solver, std::string preconditioner)
    : FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >(
            gridView, problem,
            *(new LevelCommunicate<typename GridView::Grid>(gridView.grid(), gridView.grid().maxLevel())),
            gridView.grid().maxLevel(), solver, preconditioner)
    {}
};

template<class GridView, class Scalar, class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = GroundwaterEquationLocalStiffness<GridView, Scalar, Problem> >
class LevelFEPressure2P: public FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid>  >
{
public:
    LevelFEPressure2P(GridView& gridView, Problem& problem, int level = -1)
    : FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >(
            gridView,
            problem,
            *(new LevelCommunicate<typename GridView::Grid>(gridView.grid(), level)),
            level)
    {}
};

template<class GridView, class Scalar, class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = GroundwaterEquationLocalStiffness<GridView, Scalar, Problem> >
class LeafFEPressure2P: public FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LeafCommunicate<typename GridView::Grid>  >
{
public:
    LeafFEPressure2P(GridView& gridView, Problem& problem, std::string solver = "CG", std::string preconditioner = "SeqILU0")
    : FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LeafCommunicate<typename GridView::Grid> >(
            gridView,
            problem,
            *(new LeafCommunicate<typename GridView::Grid>(gridView.grid())), -1,
            solver, preconditioner)
    {}
};


}

#include "fetotalvelocity2p.hh"

#endif
