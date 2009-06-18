// $Id: mimeticpressure2p.hh 2143 2009-06-17 18:21:10Z bernd $
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *                            *
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
template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
class MimeticPressure2PBase: public Diffusion<GridView, Scalar, VC, Problem>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        Sw = 0, Sn = 1, other = 999
    };
typedef    typename GridView::Grid Grid;
    typedef CRFunction<Grid,Scalar,GridView,Communication,1> TraceType;
    typedef P0Function<GridView,Scalar,2*GridView::dimension> NormalVelType;
    typedef MimeticOperatorAssembler<Grid,Scalar,GridView,Communication,1> OperatorAssembler;
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
        DUNE_THROW(Dune::NotImplemented, "upscaled velocities only implemented in MimeticDiffusion");
    }

    //constitutive functions are initialized and stored in the variables object
    void initializeMaterialLaws();

    void vtkout (const char* name, int k) const
    {
        this->diffProblem().variables().vtkout(name, k);
    }

    MimeticPressure2PBase(GridView& gridView, Problem& prob, Communication& comm,
    		std::string satType, int level, bool calcPressure, std::string solver, std::string preconditioner)
    : Diffusion<GridView, Scalar, VC, Problem>(gridView, prob),
    saturationType((satType == "Sw") ? Sw : ((satType == "Sn") ? Sn : other)),
    level_((level >= 0) ? level : gridView.grid().maxLevel()),
    pressTrace(gridView.grid(), gridView, comm), normalVelocity(gridView),
    f(gridView.grid(), gridView, comm), A(gridView.grid(), gridView, comm),
    solverName_(solver), preconditionerName_(preconditioner)
    {
        *pressTrace = 0;
        *f = 0;
        if (saturationType == other)
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
    OperatorAssembler A;
    std::string solverName_;
    std::string preconditionerName_;
};

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
void MimeticPressure2PBase<GridView, Scalar, VC, Problem,LocalStiffnessType,Communication>::solve()
{
	std::cout << "MimeticPressure2P: solve for pressure" << std::endl;

    typedef typename CRFunction<Grid,Scalar,GridView,Communication,1>::RepresentationType Vector;
    typedef typename CROperatorAssembler<Grid,Scalar,GridView,Communication,1>::RepresentationType Matrix;
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
    		solver.apply(*pressTrace, *f, result);
    	}
    	else if (solverName_ == "BiCGSTAB")
    	{
    		BiCGSTABSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
    		solver.apply(*pressTrace, *f, result);
    	}
    	else
    		DUNE_THROW(NotImplemented, "MimeticPressure2P :: solve : combination "
    				<< preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
    	SeqPardiso<Matrix,Vector,Vector> preconditioner(*A);
    	if (solverName_ == "Loop")
    	{
    		LoopSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
    		solver.apply(*pressTrace, *f, result);
    	}
    	else
    		DUNE_THROW(NotImplemented, "MimeticPressure2P :: solve : combination "
    				<< preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
    	DUNE_THROW(NotImplemented, "MimeticPressure2P :: solve : preconditioner "
    			<< preconditionerName_ << ".");

    return;
}

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
void MimeticPressure2PBase<GridView, Scalar, VC, Problem,LocalStiffnessType,Communication>::calculateVelocity(const Scalar t=0) const
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

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
void MimeticPressure2PBase<GridView, Scalar, VC, Problem,LocalStiffnessType,Communication>::initializeMaterialLaws()
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

template<
        class GridView,
        class Scalar,
        class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = Dune::MimeticGroundwaterEquationLocalStiffness<GridView,Scalar,VC, Problem> >
class MimeticPressure2P: public MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >
{
public:
    MimeticPressure2P(GridView& gridView, Problem& problem, std::string satType = "Sw",
    		int level = -1, bool calcPressure = true, std::string solver = "CG", std::string preconditioner = "SeqILU0")
    : MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >(
            gridView, problem,
            *(new LevelCommunicate<typename GridView::Grid>(gridView.grid(), level)),
            satType, level, calcPressure, solver, preconditioner)
    {}

    MimeticPressure2P(GridView& gridView, Problem& problem, std::string solver, std::string preconditioner)
    : MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >(
            gridView, problem,
            *(new LevelCommunicate<typename GridView::Grid>(gridView.grid(), gridView.grid().maxLevel())),
            "Sw", gridView.grid().maxLevel(), true, solver, preconditioner)
    {}
};

template<
        class GridView,
        class Scalar,
        class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = Dune::MimeticGroundwaterEquationLocalStiffness<GridView,Scalar,VC, Problem> >
class LevelMimeticPressure2P: public MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >
{
public:
    LevelMimeticPressure2P(GridView& gridView, Problem& problem, std::string satType = "Sw", int level = -1, bool calcPressure = true)
    : MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LevelCommunicate<typename GridView::Grid> >(
            gridView,
            problem,
            *(new LevelCommunicate<typename GridView::Grid>(gridView.grid(), level)),
            satType,
            level,
            calcPressure)
    {}
};

template<
        class GridView,
        class Scalar,
        class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> ,
        class LocalStiffnessType = Dune::MimeticGroundwaterEquationLocalStiffness<GridView,Scalar,VC, Problem> >
class LeafMimeticPressure2P: public MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LeafCommunicate<typename GridView::Grid> >
{
public:
    LeafMimeticPressure2P(GridView& gridView, Problem& problem, std::string satType = "Sw", int level = -1, bool calcPressure = true)
    : MimeticPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, LeafCommunicate<typename GridView::Grid> >(
            gridView,
            problem,
            *(new LeafCommunicate<typename GridView::Grid>(gridView.grid())),
            satType,
            level,
            calcPressure)
    {}
};

}
#endif
