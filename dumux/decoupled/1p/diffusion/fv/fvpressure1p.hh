// $Id: fvpressure2p.hh 3826 2010-07-14 07:03:41Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Markus Wolff                                 *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVPRESSURE1P_HH
#define DUMUX_FVPRESSURE1P_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include "dumux/common/pardiso.hh"
#include <dumux/decoupled/1p/1pproperties.hh>

/**
 * @file
 * @brief  Single Phase Finite Volume Model
 * @author Markus Wolff
 */

namespace Dumux
{

//! \ingroup OnePhase
//! \brief Single Phase Finite Volume Model
/*! Provides a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v} = q.\f]
 * The velocity \f$\boldsymbol{v}\f$ is the single phase Darcy velocity:
 * \f[ \boldsymbol{v} = -\frac{1}{\mu} \boldsymbol{K} \left(\text{grad}\, p + \rho g  \text{grad}\, z\right), \f]
 * where \f$p\f$ is the pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\mu\f$ the viscosity, \f$\rho\f$ the density, and \f$g\f$ the gravity constant,
 * and \f$q\f$ is the source term.
 * At the boundary, \f$p = p_D\f$ on \f$\Gamma_{Dirichlet}\f$, and \f$\boldsymbol{v}_{total}  = q_N\f$
 * on \f$\Gamma_{Neumann}\f$.
 *
 * @tparam TypeTag The Type Tag
 *
 */
template<class TypeTag> class FVPressure1P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Fluid)) Fluid;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    //function which assembles the system of equations to be solved
    void assemble(bool first);

    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

protected:
    //! Returns reference to the instance of the problem definition
    Problem& problem()
    {
        return problem_;
    }
    //! Returns reference to the instance of the problem definition
    const Problem& problem() const
    {
        return problem_;
    }

public:
    //! Initializes the problem
    /*!
     *  @param solveTwice repeats the pressure calculation step
     *
     *  Calculates the pressure \f$p\f$ as solution of the boundary value
     *  \f[  \text{div}\, \boldsymbol{v} = q, \f]
     *  subject to appropriate boundary conditions.
     */

    void initialize(bool solveTwice = true)
    {
        assemble(true);
        solve();
        if (solveTwice)
        {
            assemble(false);
            solve();
        }
        return;
    }

    //! Calculates the pressure.
    /*!
     *  @param solveTwice without any function here!
     *
     *  Calculates the pressure \f$p\f$ as solution of the boundary value
     *  \f[  \text{div}\, \boldsymbol{v} = q, \f]
     *  subject to appropriate boundary conditions.
     */
    void pressure(bool solveTwice = true)
    {
        assemble(false);
        solve();

        return;
    }

    // serialization methods
    //! Function needed for restart option.
    template<class Restarter>
    void serialize(Restarter &res)
    {
        return;
    }

    //! Function needed for restart option.
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        return;
    }

    //! \brief Writes data files
    /*  \param writer VTK-Writer for the current simulation run */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        typename Variables::ScalarSolutionType *pressure = writer.template createField<Scalar, 1> (
                problem_.gridView().size(0));

        *pressure = problem_.variables().pressure();

        writer.addCellData(pressure, "pressure");

        return;
    }

    //! Constructs a FVPressure1P object
    /**
     * \param problem a problem class object
     */
    FVPressure1P(Problem& problem) :
        problem_(problem), A_(problem.variables().gridSize(), problem.variables().gridSize(), (2 * dim + 1)
                * problem.variables().gridSize(), Matrix::random), f_(problem.variables().gridSize()), gravity(
                problem.gravity())
    {
        initializeMatrix();
    }

private:
    Problem& problem_;
    Matrix A_;
    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > f_;
protected:
    const Dune::FieldVector<Scalar, dimWorld>& gravity; //!< vector including the gravity constant
};

//!initializes the matrix to store the system of equations
template<class TypeTag>
void FVPressure1P<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
                rowSize++;
        A_.setrowsize(globalIdxI, rowSize);
    }
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem_.variables().index(*outside);

                // add off diagonal index
                A_.addindex(globalIdxI, globalIdxJ);
            }
    }
    A_.endindices();

    return;
}

//!function which assembles the system of equations to be solved
template<class TypeTag>
void FVPressure1P<TypeTag>::assemble(bool first)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_ = 0;

    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().center();

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        Scalar temperatureI = problem_.temperature(globalPos, *eIt);
        Scalar referencePressI = problem_.referencePressure(globalPos, *eIt);

        Scalar densityI = Fluid::density(temperatureI, referencePressI);

        // set right side to zero
        Scalar source = problem_.source(globalPos, *eIt) / densityI;

        f_[globalIdxI] = volume * source;

        // get absolute permeability
        FieldMatrix permeabilityI(problem_.spatialParameters().intrinsicPermeability(globalPos, *eIt));

        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            int isIndex = isIt->indexInInside();

            // get normal vector
            Dune::FieldVector < Scalar, dimWorld > unitOuterNormal = isIt->centerUnitOuterNormal();

            // get face volume
            Scalar faceArea = isIt->geometry().volume();

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector < Scalar, dimWorld > distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                FieldMatrix permeabilityJ = problem_.spatialParameters().intrinsicPermeability(globalPosNeighbor,
                        *neighborPointer);

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);

                // harmonic mean of permeability
                for (int x = 0; x < dim; x++)
                {
                    meanPermeability[x][x] = 2 * permeabilityI[x][x] * permeabilityJ[x][x] / (permeabilityI[x][x]
                            + permeabilityJ[x][x]);
                    for (int y = 0; y < dim; y++)
                    {
                        if (x != y)
                        {//use arithmetic mean for the off-diagonal entries to keep the tensor property!
                            meanPermeability[x][y] = 0.5 * (permeabilityI[x][y] + permeabilityJ[x][y]);
                        }
                    }
                }

                Dune::FieldVector < Scalar, dim > permeability(0);
                meanPermeability.mv(unitOuterNormal, permeability);

                Scalar temperatureJ = problem_.temperature(globalPosNeighbor, *neighborPointer);
                Scalar referencePressJ = problem_.referencePressure(globalPosNeighbor, *neighborPointer);

                Scalar densityJ = Fluid::density(temperatureJ, referencePressJ);

                Scalar rhoMean = 0.5 * (densityI + densityJ);

                // update diagonal entry
                Scalar entry;

                //calculate potential gradients
                Scalar potential = 0;

                Scalar density = 0;

                //if we are at the very first iteration we can't calculate phase potentials
                if (!first)
                {
                    potential = problem_.variables().potential(globalIdxI, isIndex);

                    density = (potential > 0.) ? densityI : densityJ;

                    density = (potential == 0.) ? rhoMean : density;

                    potential = (problem_.variables().pressure()[globalIdxI]
                            - problem_.variables().pressure()[globalIdxJ]) / dist;

                    potential += density * (unitOuterNormal * gravity);

                    //store potentials for further calculations (velocity, saturation, ...)
                    problem_.variables().potential(globalIdxI, isIndex) = potential;
                }

                //do the upwinding depending on the potentials

                density = (potential > 0.) ? densityI : densityJ;

                density = (potential == 0) ? rhoMean : density;

                //calculate current matrix entry
                entry = ((permeability * unitOuterNormal) / dist) * faceArea;

                //calculate right hand side
                Scalar rightEntry = density * (permeability * gravity) * faceArea;

                //set right hand side
                f_[globalIdxI] -= rightEntry;

                // set diagonal entry
                A_[globalIdxI][globalIdxI] += entry;

                // set off-diagonal entry
                A_[globalIdxI][globalIdxJ] = -entry;
            }

            // boundary face

            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                Dune::FieldVector < Scalar, dimWorld > distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = problem_.bctype(globalPosFace, *isIt);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    //permeability vector at boundary
                    Dune::FieldVector < Scalar, dim > permeability(0);
                    permeabilityI.mv(unitOuterNormal, permeability);

                    //get dirichlet pressure boundary condition
                    Scalar pressBound = problem_.dirichlet(globalPosFace, *isIt);

                    //calculate current matrix entry
                    Scalar entry = ((permeability * unitOuterNormal) / dist) * faceArea;

                    //calculate right hand side
                    Scalar rightEntry = densityI * (permeability * gravity) * faceArea;

                    // set diagonal entry and right hand side entry
                    A_[globalIdxI][globalIdxI] += entry;
                    f_[globalIdxI] += entry * pressBound;
                    f_[globalIdxI] -= rightEntry;
                }
                //set neumann boundary condition

                else
                {
                    Scalar J = problem_.neumann(globalPosFace, *isIt) / densityI;

                    f_[globalIdxI] -= J * faceArea;
                }
            }
        } // end all intersections

    } // end grid traversal
    return;
}

//!solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure1P<TypeTag>::solve()
{
    typedef typename GET_PROP(TypeTag, PTAG(SolverParameters)) SolverParameters;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressurePreconditioner)) Preconditioner;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureSolver)) Solver;

    Dune::MatrixAdapter < Matrix, Vector, Vector > op(A_); // make linear operator from A_
    Dune::InverseOperatorResult result;

    double reduction = SolverParameters::reductionSolver;
    int maxItSolver = SolverParameters::maxIterationNumberSolver;
    int iterPreconditioner = SolverParameters::iterationNumberPreconditioner;
    int verboseLevelSolver = SolverParameters::verboseLevelSolver;
    double relaxation = SolverParameters::relaxationPreconditioner;

    if (verboseLevelSolver)
        std::cout << "FVPressure1P: solve for pressure" << std::endl;

    Preconditioner preconditioner(A_, iterPreconditioner, relaxation);
    Solver solver(op, preconditioner, reduction, maxItSolver, verboseLevelSolver);
    solver.apply(problem_.variables().pressure(), f_, result);

    //                printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
    //                printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
    //                printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);

    return;
}

}
#endif
