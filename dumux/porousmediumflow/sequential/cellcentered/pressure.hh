// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVPRESSURE_HH
#define DUMUX_FVPRESSURE_HH

// dumux environment
#include <type_traits>
#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <map>
/**
 * @file
 * @brief  Finite Volume Diffusion Model
 */

namespace Dumux
{
//! The finite volume base class for the solution of a pressure equation
/*! \ingroup IMPET
 *  Base class for finite volume (FV) implementations of a diffusion-like pressure equation.
 *  The class provides a methods for assembling of the global matrix and right hand side (RHS)
 *  as well as for solving the system of equations.
 *  Additionally, it contains the global matrix, the RHS-vector as well as the solution vector.
 *  A certain pressure equation defined in the implementation of this base class must be splitted
 *  into a storage term, a flux term and a source term.
 *  Corresponding functions (<tt>getSource()</tt>, <tt>getStorage()</tt>, <tt>getFlux()</tt> and
 *  <tt>getFluxOnBoundary()</tt>) have to be defined in the implementation.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure
{
    //the model implementation
    using Implementation = GetPropType<TypeTag, Properties::PressureModel>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    // using declarations to abbreviate several dune classes...
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    // the typenames used for the stiffness matrix and solution vector
    using Matrix = GetPropType<TypeTag, Properties::PressureCoefficientMatrix>;
    using RHSVector = GetPropType<TypeTag, Properties::PressureRHSVector>;
    using PressureSolution = GetPropType<TypeTag, Properties::PressureSolutionVector>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

protected:

    /*! Type of the vector of entries
     *
     * Contains the return values of the get*()-functions (matrix or right-hand side entry).
     */
    using EntryType = Dune::FieldVector<Scalar, 2>;

    //! Indices of matrix and rhs entries
    /**
    * During the assembling of the global system of equations get-functions are called (getSource(),
    * getFlux(), etc.), which return global matrix or right hand side entries in a vector.
    * These can be accessed using following indices:
    */
    enum
    {
        rhs = 1,//!<index for the right hand side entry
        matrix = 0//!<index for the global matrix entry

    };

    enum
    {
        pressEqIdx = Indices::pressureEqIdx,
    };

    //!Initialize the global matrix of the system of equations to solve
    void initializeMatrix();
    void initializeMatrixRowSize();
    void initializeMatrixIndices();


    /*!\brief Function which assembles the system of equations to be solved
     *
     *  This function assembles the Matrix and the right hand side (RHS) vector to solve for
     * a pressure field with a Finite-Volume (FV) discretization.
     * Implementations of this base class have to provide the methods <tt>getSource()</tt>,
     * <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt> if the assemble() method is called!
     *
     * \param first Indicates if function is called at the initialization step or during the simulation
     * (If <tt>first</tt> is <tt>true</tt>, no pressure field of previous iterations is required)
     */
    void assemble(bool first);

    //!Solves the global system of equations to get the spatial distribution of the pressure
    void solve();

    //!Returns the vector containing the pressure solution
    PressureSolution& pressure()
    {   return pressure_;}

    //!Returns the vector containing the pressure solution
    const PressureSolution& pressure() const
    {   return pressure_;}

    //!Initialization of the pressure solution vector: Initialization with meaningful values may
    //result in better convergence of the linear solver!
    void initializePressure()
    {
        for (const auto& element : elements(problem_.gridView()))
        {
            PrimaryVariables initValues;
            problem_.initial(initValues, element);

            int eIdxGlobal = problem_.variables().index(element);

            pressure_[eIdxGlobal] = initValues[pressEqIdx];
        }
    }

public:
    /*! \brief Function which calculates the source entry
     *
     * Function computes the source term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param element Grid element
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getSource(EntryType& entry, const Element& element, const CellData& cellData, const bool first);

    /*! \brief Function which calculates the storage entry
     *
     * Function computes the storage term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param element Grid element
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getStorage(EntryType& entry, const Element& element, const CellData& cellData, const bool first);

    /*! \brief Function which calculates the flux entry
     *
     * Function computes the inter-cell flux term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param intersection Intersection of two grid elements
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getFlux(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first);

    /*! \brief Function which calculates the boundary flux entry
     *
     * Function computes the boundary-flux term and writes it to the corresponding entry of the entry vector
     *
     * \param entry Vector containing return values of the function
     * \param intersection Intersection of two grid elements
     * \param cellData Object containing all model relevant cell data
     * \param first Indicates if function is called in the initialization step or during the simulation
     */
    void getFluxOnBoundary(EntryType& entry,
            const Intersection& intersection, const CellData& cellData, const bool first);

    /*! \brief Public access function for the primary pressure variable
     *
     * Function returns the cell pressure value at index <tt>eIdxGlobal</tt>
     *
     * \param eIdxGlobal Global index of a grid cell
     */
    const Scalar pressure(const int eIdxGlobal) const
    {   return pressure_[eIdxGlobal];}

    //!Returns the global matrix of the last pressure solution step.
    const Matrix& globalMatrix()
    {
        return A_;
    }

    //!Returns the right hand side vector of the last pressure solution step.
    const RHSVector& rightHandSide()
    {
        return f_;
    }

    /*! \brief Initialize pressure model
     *
     * Function initializes the sparse matrix to solve the global system of equations and sets/calculates the initial pressure
     */
    void initialize()
    {
        int size = problem_.gridView().size(0);//resize to make sure the final grid size (after the problem was completely built) is used!
        A_.setSize(size, size);
        A_.setBuildMode(Matrix::random);
        f_.resize(size);
        pressure_.resize(size);
        initializePressure();
        asImp_().initializeMatrix();// initialize sparse matrix
    }

    /*! \brief Pressure update
     *
     * Function reassembles the system of equations and solves for a new pressure solution.
     */
    void update()
    {
        asImp_().assemble(false); Dune::dinfo << "pressure calculation"<< std::endl;
        solve();

        return;
    }

    void calculateVelocity()
    {
        DUNE_THROW(Dune::NotImplemented,"Velocity calculation not implemented in pressure model!");
    }

    void updateVelocity()
    {} //empty function for the case the velocity is calculated in the transport model

    /*! \brief  Function for serialization of the pressure field.
     *
     *  Function needed for restart option. Writes the pressure of a grid element to a restart file.
     *
     *  \param outstream Stream into the restart file.
     *  \param element Grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        outstream << pressure_[eIdxGlobal][0];
    }

    /*! \brief  Function for deserialization of the pressure field.
     *
     *  Function needed for restart option. Reads the pressure of a grid element from a restart file.
     *
     *  \param instream Stream from the restart file.
     *  \param element Grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        instream >> pressure_[eIdxGlobal][0];
    }

    /*! \brief Set a pressure to be fixed at a certain cell.
     *
     *Allows to fix a pressure somewhere (at one certain cell) in the domain.
     *This can be necessary e.g. if only Neumann boundary conditions are defined.
     *The pressure is fixed until the <tt>unsetFixPressureAtIndex()</tt> function is called
     *
     * \param pressure Pressure value at eIdxGlobal
     * \param eIdxGlobal Global index of a grid cell
     */
    void setFixPressureAtIndex(Scalar pressure, int eIdxGlobal)
    {
        fixPressure_.insert(std::make_pair(eIdxGlobal, pressure));
    }

    /*! \brief Reset the fixed pressure state
     *
     * No pressure is fixed inside the domain until <tt>setFixPressureAtIndex()</tt> function is called again.
     *
     * \param eIdxGlobal Global index of a grid cell
     */
    void unsetFixPressureAtIndex(int eIdxGlobal)
    {
        fixPressure_.erase(eIdxGlobal);
    }

    void resetFixPressureAtIndex()
    {
        fixPressure_.clear();
    }

    /*! \brief Constructs a FVPressure object
     *
     * \param problem A problem class object
     */
    FVPressure(Problem& problem) :
    problem_(problem)
    {}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}

    Problem& problem_;

    PressureSolution pressure_;

    std::string solverName_;
    std::string preconditionerName_;
protected:
    Matrix A_;//!<Global stiffnes matrix (sparse matrix which is build by the <tt> initializeMatrix()</tt> function)
    RHSVector f_;//!<Right hand side vector
private:
    std::map<int, Scalar> fixPressure_;
};

//!Initialize the global matrix of the system of equations to solve
template<class TypeTag>
void FVPressure<TypeTag>::initializeMatrix()
{
    initializeMatrixRowSize();
    A_.endrowsizes();
    initializeMatrixIndices();
    A_.endindices();
}

//!Initialize the global matrix of the system of equations to solve
template<class TypeTag>
void FVPressure<TypeTag>::initializeMatrixRowSize()
{
    // determine matrix row sizes
    for (const auto& element : elements(problem_.gridView()))
    {
        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        for (const auto& intersection : intersections(problem_.gridView(), element))
        {
            if (intersection.neighbor())
                rowSize++;
        }
        A_.setrowsize(eIdxGlobalI, rowSize);
    }
}

//!Initialize the global matrix of the system of equations to solve
template<class TypeTag>
void FVPressure<TypeTag>::initializeMatrixIndices()
{
    // determine position of matrix entries
    for (const auto& element : elements(problem_.gridView()))
    {
        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        // add diagonal index
        A_.addindex(eIdxGlobalI, eIdxGlobalI);

        // run through all intersections with neighbors
        for (const auto& intersection : intersections(problem_.gridView(), element))
            if (intersection.neighbor())
            {
                // access neighbor
                int eIdxGlobalJ = problem_.variables().index(intersection.outside());

                // add off diagonal index
                A_.addindex(eIdxGlobalI, eIdxGlobalJ);
            }
    }
}


/*!\brief Function which assembles the system of equations to be solved
 *
 *  This function assembles the Matrix and the right hand side (RHS) vector to solve for
 * a pressure field with a Finite-Volume (FV) discretization.
 * Implementations of this base class have to provide the methods <tt>getSource()</tt>,
 * <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt> if the assemble() method is called!
 *
 * \param first Indicates if function is called at the initialization step or during the simulation
 *              (If <tt>first</tt> is <tt>true</tt>, no pressure field of previous iterations is required)
 */
template<class TypeTag>
void FVPressure<TypeTag>::assemble(bool first)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_ = 0;

    for (const auto& element : elements(problem_.gridView()))
    {
        // get the global index of the cell
        int eIdxGlobalI = problem_.variables().index(element);

        // assemble interior element contributions
        if (element.partitionType() == Dune::InteriorEntity)
        {
            // get the cell data
            CellData& cellDataI = problem_.variables().cellData(eIdxGlobalI);

            EntryType entries(0.);

            /*****  source term ***********/
            asImp_().getSource(entries, element, cellDataI, first);
            f_[eIdxGlobalI] += entries[rhs];

            /*****  flux term ***********/
            // iterate over all faces of the cell
            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                /************* handle interior face *****************/
                if (intersection.neighbor())
                {
                    auto elementNeighbor = intersection.outside();

                    int eIdxGlobalJ = problem_.variables().index(elementNeighbor);

                    // check for hanging nodes
                    bool haveSameLevel = (element.level() == elementNeighbor.level());
                    // calculate only from one side (except for hanging nodes), but add matrix entries for both sides
                    // the last condition is needed to properly assemble in the presence
                    // of ghost elements
                    if (getPropValue<TypeTag, Properties::VisitFacesOnlyOnce>()
                        && (eIdxGlobalI > eIdxGlobalJ) && haveSameLevel
                        && elementNeighbor.partitionType() == Dune::InteriorEntity)
                        continue;

                    entries = 0;
                    asImp_().getFlux(entries, intersection, cellDataI, first);

                    //set right hand side
                    f_[eIdxGlobalI] -= entries[rhs];

                    // set diagonal entry
                    A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];

                    // set off-diagonal entry
                    A_[eIdxGlobalI][eIdxGlobalJ] -= entries[matrix];

                    // The second condition is needed to not spoil the ghost element entries
                    if (getPropValue<TypeTag, Properties::VisitFacesOnlyOnce>()
                        && elementNeighbor.partitionType() == Dune::InteriorEntity)
                    {
                        f_[eIdxGlobalJ] += entries[rhs];
                        A_[eIdxGlobalJ][eIdxGlobalJ] += entries[matrix];
                        A_[eIdxGlobalJ][eIdxGlobalI] -= entries[matrix];
                    }

                } // end neighbor

                /************* boundary face ************************/
                else
                {
                    entries = 0;
                    asImp_().getFluxOnBoundary(entries, intersection, cellDataI, first);

                    //set right hand side
                    f_[eIdxGlobalI] += entries[rhs];
                    // set diagonal entry
                    A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
                }
            } //end interfaces loop
    //        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);

            /*****  storage term ***********/
            entries = 0;
            asImp_().getStorage(entries, element, cellDataI, first);
            f_[eIdxGlobalI] += entries[rhs];
    //         set diagonal entry
            A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
        }
        // assemble overlap and ghost element contributions
        else
        {
            A_[eIdxGlobalI] = 0.0;
            A_[eIdxGlobalI][eIdxGlobalI] = 1.0;
            f_[eIdxGlobalI] = pressure_[eIdxGlobalI];
        }
    } // end grid traversal
//    printmatrix(std::cout, A_, "global stiffness matrix after assempling", "row", 11,3);
//    printvector(std::cout, f_, "right hand side", "row", 10);
}

// forward declaration
template<class T>
class AMGBiCGSTABBackend;

namespace Detail {
template<class T> struct isParallelAMGBackend : public std::false_type {};
template<class T>
struct isParallelAMGBackend<Dumux::AMGBiCGSTABBackend<T>> : public std::true_type {};
} // end namespace Detail

template<class Solver, class Problem>
typename std::enable_if_t<!Detail::isParallelAMGBackend<Solver>::value, Solver>
getSolver(const Problem& problem)
{
    return Solver();
}

template<class Solver, class Problem>
typename std::enable_if_t<Detail::isParallelAMGBackend<Solver>::value, Solver>
getSolver(const Problem& problem)
{
    return Solver(problem.gridView(), problem.model().dofMapper());
}

//!Solves the global system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure<TypeTag>::solve()
{
    using Solver = GetPropType<TypeTag, Properties::LinearSolver>;

    int verboseLevelSolver = getParam<int>("LinearSolver.Verbosity", 0);

    if (verboseLevelSolver)
        std::cout << __FILE__ << ": solve for pressure" << std::endl;

    //set a fixed pressure for a certain cell
    if (fixPressure_.size() > 0)
    {
        for (auto it = fixPressure_.begin(); it != fixPressure_.end(); ++it)
        {
            A_[it->first] = 0;
            A_[it->first][it->first] = 1;
            f_[it->first] = it->second;
        }
    }

//    printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//    printvector(std::cout, f_, "right hand side", "row", 10, 1, 3);

    auto solver = getSolver<Solver>(problem_);
    bool converged = solver.solve(A_, pressure_, f_);

    if (!converged)
        DUNE_THROW(Dumux::NumericalProblem, "Pressure solver did not converge!");

//    printvector(std::cout, pressure_, "pressure", "row", 200, 1, 3);
}

} //end namespace Dumux
#endif
