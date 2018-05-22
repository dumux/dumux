// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
#ifndef DUMUX_IMPLICIT_ASSEMBLER_HH
#define DUMUX_IMPLICIT_ASSEMBLER_HH

#include "properties.hh"

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class ImplicitAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum{ dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

public:

    //! copying the jacobian assembler is not a good idea
    ImplicitAssembler(const ImplicitAssembler&) = delete;

    //! default constructor
    ImplicitAssembler() : problemPtr_(nullptr)
    {}

    /*!
     * \brief Initialize the jacobian assembler.
     *
     * At this point we can assume that all objects in the problem and
     * the model have been allocated. We can not assume that they are
     * fully initialized, though.
     *
     * \param problem The problem object
     */
    void init(Problem& problem)
    {
        problemPtr_ = &problem;

        // initialize the BCRS matrix
        asImp_().createMatrix_();

        // initialize the jacobian matrix with zeros
        *matrix_ = 0;

        // allocate the residual vector
        residual_.resize(problem.model().numDofs());
    }

    /*!
     * \brief Assemble the global Jacobian of the residual and the residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        bool succeeded;
        try {
            asImp_().assemble_();
            succeeded = true;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << problem_().gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            DUNE_THROW(NumericalProblem,
                       "A process did not succeed in linearizing the system");
        }
        // printmatrix(std::cout, matrix(), "", "");

    }

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     */
    const JacobianMatrix& matrix() const
    { return *matrix_; }
    JacobianMatrix& matrix()
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const SolutionVector& residual() const
    { return residual_; }
    SolutionVector& residual()
    { return residual_; }

protected:

    // reset the global linear system of equations. if partial
    // reassemble is enabled, this means that the jacobian matrix must
    // only be erased partially!
    void resetSystem_()
    {
        // reset the right hand side.
        residual_ = 0.0;

        // reset the matrix
        (*matrix_) = 0;
    }

    // linearize the whole system
    void assemble_()
    {
        asImp_().resetSystem_();

//	std::cout << std::endl;
//	printmatrix(std::cout, this->matrix(), "Matrix", "", 7, 0);
        // assemble the elements...
        for (const auto& element : elements(gridView_()))
            this->model_().localJacobian().assemble(element, this->matrix(), this->residual());
    }

protected:
    Problem &problem_()
    { return *problemPtr_; }
    const Problem &problem_() const
    { return *problemPtr_; }
    const Model &model_() const
    { return problem_().model(); }
    Model &model_()
    { return problem_().model(); }
    const GridView &gridView_() const
    { return problem_().gridView(); }
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); }
    const ElementMapper &elementMapper_() const
    { return problem_().elementMapper(); }

    Problem *problemPtr_;

    // the jacobian matrix
    std::shared_ptr<JacobianMatrix> matrix_;
    // the right-hand side
    SolutionVector residual_;

private:

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        auto numDofs = problem_().model().numDofs();

        // allocate raw matrix
        matrix_ = std::make_shared<JacobianMatrix>(numDofs, numDofs, JacobianMatrix::random);

        // set the row sizes
        asImp_().setRowSizes_();

        // set the indices
        asImp_().addIndices_();
    }

    //! Set the row sizes
    void setRowSizes_()
    {
        DUNE_THROW(Dune::NotImplemented, "Actual implementation does not provide a setRowSizes_() method!");
    }

    void addIndices_()
    {
        DUNE_THROW(Dune::NotImplemented, "Actual implementation does not provide a addIndices_() method!");
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
