// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_FVPRESSURE_HH
#define DUMUX_FVPRESSURE_HH

// dumux environment
#include "dumux/common/math.hh"
#include <dumux/decoupled/common/pressureproperties.hh>
/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Benjamin Faigle, Bernd Flemisch, Markus Wolff
 */

namespace Dumux
{
//! The finite volume base class for the solution of a pressure equation
/*! \ingroup IMPET
 *  Base class for finite volume (FV) implementations of a diffusion-like pressure equation.
 *  The class provides a methods for assembling of the global matrix and right hand side (RHS) as well as for solving the system of equations.
 *  Additionally, it contains the global matrix, the RHS-vector as well as the solution vector.
 *  A certain pressure equation defined in the implementation of this base class must be splitted into a storage term, a flux term and a source term.
 *  Corresponding functions (<tt>getSource()</tt>, <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt>) have to be defined in the implementation.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure
{
    //the model implementation
    typedef typename GET_PROP_TYPE(TypeTag, PressureModel) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP(TypeTag, SolutionTypes)::ScalarSolution ScalarSolution;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    // typedefs to abbreviate several dune classes...
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    // the typenames used for the stiffness matrix and solution vector
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) RHSVector;

protected:

    //! Indices of matrix and rhs entries
    /**
    * During the assembling of the global system of equations get-functions are called (getSource(), getFlux(), etc.), which return global matrix or right hand side entries in a vector. These can be accessed using following indices:
    */
    enum
    {
        rhs = 1,//!<index for the right hand side entry
        matrix = 0//!<index for the global matrix entry

    };

    //!Initialize the global matrix of the system of equations to solve
    void initializeMatrix();

    //! Function which assembles the system of equations to be solved
    /* This function assembles the Matrix and the right hand side (RHS) vector to solve for
     * a pressure field with a Finite-Volume (FV) discretization. Implementations of this base class have to provide the methods
     * <tt>getSource()</tt>, <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt> if the assemble() method is called!
     *
     * \param first Indicates if function is called at the initialization step or during the simulation (If <tt>first</tt> is <tt>true</tt>, no pressure field of previous iterations is required)
     */
    void assemble(bool first);

    //!Solves the global system of equations to get the spatial distribution of the pressure
    void solve();

    //!Returns the vector containing the pressure solution
    ScalarSolution& pressure()
    {   return pressure_;}

    //!Returns the vector containing the pressure solution
    const ScalarSolution& pressure() const
    {   return pressure_;}

public:
    //! Function which calculates the source entry
    /**
     * Function computes the source term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param element Grid element
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getSource(Dune::FieldVector<Scalar, 2>& entry, const Element& element, const CellData& cellData, const bool first);

    //!Function which calculates the storage entry
    /**
     * Function computes the storage term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param element Grid element
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getStorage(Dune::FieldVector<Scalar, 2>& entry, const Element& element, const CellData& cellData, const bool first);

    //!Function which calculates the flux entry
    /**
     * Function computes the inter-cell flux term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param intersection Intersection of two grid elements
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getFlux(Dune::FieldVector<Scalar, 2>& entry, const Intersection& intersection, const CellData& cellData, const bool first);

    //!Function which calculates the boundary flux entry
    /**
     * Function computes the boundary-flux term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param intersection Intersection of two grid elements
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entry,
            const Intersection& intersection, const CellData& cellData, const bool first);

    //! Public access function for the primary pressure variable
    /**
     * Function returns the cell pressure value at index <tt>globalIdx</tt>
     *
     * \param globalIdx Global index of a grid cell
     */
    const Scalar pressure(int globalIdx) const
    {   return pressure_[globalIdx];}

    //!Initialize pressure model
    /**
     * Function initializes the sparse matrix to solve the global system of equations and sets/calculates the initial pressure
     */
    void initialize()
    {
        size_ = problem_.gridView().size(0);//resize to make sure the final grid size (after the problem was completely built) is used!
        A_.setSize(size_, size_);
        f_.resize(size_);
        pressure_.resize(size_);
        initializeMatrix();// initialize sparse matrix
        pressure_ = 0;
    }

    //!Pressure update
    /**
     * Function reassembles the system of equations and solves for a new pressure solution.
     */
    void update()
    {
        assemble(false); Dune::dinfo << "pressure calculation"<< std::endl;
        solve();

        return;
    }

    //! Function for serialization of the pressure field.
    /** Function needed for restart option. Writes the pressure of a grid element to a restart file.
     *
     *  \param outstream Stream into the restart file.
     *  \param element Grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int globalIdx = problem_.variables().index(element);
        outstream << pressure_[globalIdx][0];
    }

    //! Function for deserialization of the pressure field.
    /** Function needed for restart option. Reads the pressure of a grid element from a restart file.
     *
     *  \param instream Stream from the restart file.
     *  \param element Grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = problem_.variables().index(element);
        instream >> pressure_[globalIdx][0];
    }

    //! Set a pressure to be fixed at a certain cell.
    /**
     *Allows to fix a pressure somewhere (at one certain cell) in the domain. This can be necessary e.g. if only Neumann boundary conditions are defined.
     *The pressure is fixed until the <tt>unsetFixPressureAtIndex()</tt> function is called
     *
     * \param pressure Pressure value at globalIdx
     * \param globalIdx Global index of a grid cell
     */
    void setFixPressureAtIndex(Scalar pressure, int globalIdx)
    {
        hasFixPressureAtIndex_ = true;
        fixPressureAtIndex_ = pressure;
        idxFixPressureAtIndex_ = globalIdx;
    }

    //! Reset the fixed pressure state
    /**
     * No pressure is fixed inside the domain until <tt>setFixPressureAtIndex()</tt> function is called again.
     *
     * \param globalIdx Global index of a grid cell
     */
    void unsetFixPressureAtIndex(int globalIdx)
    {
        hasFixPressureAtIndex_ = false;
        fixPressureAtIndex_ = 0.0;
        idxFixPressureAtIndex_ = 0.0;
    }

    //! Constructs a FVPressure object
    /**
     * \param problem A problem class object
     */
    FVPressure(Problem& problem) :
    problem_(problem), size_(problem.gridView().size(0)),
    pressure_(size_), A_(size_, size_, Matrix::random), f_(size_),
    fixPressureAtIndex_(0),
    idxFixPressureAtIndex_(0),
    hasFixPressureAtIndex_(false)
    {}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}

    Problem& problem_;

    int size_;
    ScalarSolution pressure_;

    std::string solverName_;
    std::string preconditionerName_;
protected:
    Matrix A_;//!<Global stiffnes matrix (sparse matrix which is build by the <tt> initializeMatrix()</tt> function)
    RHSVector f_;//!<Right hand side vector
private:
    Scalar fixPressureAtIndex_;
    Scalar idxFixPressureAtIndex_;
    bool hasFixPressureAtIndex_;
};

//!Initialize the global matrix of the system of equations to solve
template<class TypeTag>
void FVPressure<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            if (isIt->neighbor())
                rowSize++;
        }
        A_.setrowsize(globalIdxI, rowSize);
    }
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
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

//! Function which assembles the system of equations to be solved
/* This function assembles the Matrix and the right hand side (RHS) vector to solve for
 * a pressure field with a Finite-Volume (FV) discretization. Implementations of this base class have to provide the methods
 * <tt>getSource()</tt>, <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt> if the assemble() method is called!
 *
 * \param first Indicates if function is called at the initialization step or during the simulation (If <tt>first</tt> is <tt>true</tt>, no pressure field of previous iterations is required)
 */
template<class TypeTag>
void FVPressure<TypeTag>::assemble(bool first)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_ = 0;

    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell information
        int globalIdxI = problem_.variables().index(*eIt);
        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        Dune::FieldVector<Scalar, 2> entries(0.);

        /*****  source term ***********/
        asImp_().getSource(entries, *eIt, cellDataI, first);
        f_[globalIdxI] += entries[rhs];

        /*****  flux term ***********/
        // iterate over all faces of the cell
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            /************* handle interior face *****************/
            if (isIt->neighbor())
            {
                ElementPointer elementNeighbor = isIt->outside();

                int globalIdxJ = problem_.variables().index(*elementNeighbor);

                //check for hanging nodes
                //take a hanging node never from the element with smaller level!
                bool haveSameLevel = (eIt->level() == elementNeighbor->level());
                //calculate only from one side, but add matrix entries for both sides
                if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce) && (globalIdxI > globalIdxJ) && haveSameLevel)
                    continue;

                //check for hanging nodes
                entries = 0;
                asImp_().getFlux(entries, *isIt, cellDataI, first);

                //set right hand side
                f_[globalIdxI] -= entries[rhs];

                // set diagonal entry
                A_[globalIdxI][globalIdxI] += entries[matrix];

                    // set off-diagonal entry
                A_[globalIdxI][globalIdxJ] -= entries[matrix];

                if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce))
                {
                    f_[globalIdxJ] += entries[rhs];
                    A_[globalIdxJ][globalIdxJ] += entries[matrix];
                    A_[globalIdxJ][globalIdxI] -= entries[matrix];
                }

            } // end neighbor

            /************* boundary face ************************/
            else
            {
                entries = 0;
                asImp_().getFluxOnBoundary(entries, *isIt, cellDataI, first);

                //set right hand side
                f_[globalIdxI] += entries[rhs];
                // set diagonal entry
                A_[globalIdxI][globalIdxI] += entries[matrix];
            }
        } //end interfaces loop
//        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);

        /*****  storage term ***********/
        entries = 0;
        asImp_().getStorage(entries, *eIt, cellDataI, first);
        f_[globalIdxI] += entries[rhs];
//         set diagonal entry
        A_[globalIdxI][globalIdxI] += entries[matrix];
    } // end grid traversal
//    printmatrix(std::cout, A_, "global stiffness matrix after assempling", "row", 11,3);
//    printvector(std::cout, f_, "right hand side", "row", 10);
    return;
}

//!Solves the global system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure<TypeTag>::solve()
{
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) Solver;

    int verboseLevelSolver = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);

    if (verboseLevelSolver)
        std::cout << __FILE__ << ": solve for pressure" << std::endl;

    //set a fixed pressure for a certain cell
    if (hasFixPressureAtIndex_)
    {
        A_[idxFixPressureAtIndex_] = 0;
        A_[idxFixPressureAtIndex_][idxFixPressureAtIndex_] = 1;
        f_[idxFixPressureAtIndex_] = fixPressureAtIndex_;
    }

//    printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//    printvector(std::cout, f_, "right hand side", "row", 10, 1, 3);

    Solver solver(problem_);
    solver.solve(A_, pressure_, f_);

//    printvector(std::cout, pressure_, "pressure", "row", 200, 1, 3);

    return;
}

} //end namespace Dumux
#endif
