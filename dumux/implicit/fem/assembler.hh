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
#ifndef DUMUX_FEM_ASSEMBLER_HH
#define DUMUX_FEM_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/assembler.hh>

#include <dune/istl/io.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief A finite element assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class FemAssembler : public ImplicitAssembler<TypeTag>
{
    using ParentType = ImplicitAssembler<TypeTag>;
    friend ParentType;

    using FEBasis = typename GET_PROP_TYPE(TypeTag, FeBasis);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);

    using ctype = typename GridView::ctype;
    using ReferenceElements = typename Dune::ReferenceElements<ctype, dim>;

public:
    void init(Problem& problem)
    {
        feBasisPtr_ = std::make_shared<FEBasis>(problem.gridView());
        ParentType::init(problem);
    }

    const FEBasis& feBasis() const
    { return *feBasisPtr_; }

    std::size_t numDofs() const
    { return feBasis().size(); }

protected:

    // linearize the whole system
    void assemble_()
    {
//    	std::cout << std::endl;
//    	printmatrix(std::cout, this->matrix(), "vorAssembleMatrix", "", 7, 0);
        ParentType::assemble_();

//        std::cout << std::endl;
//        printmatrix(std::cout, this->matrix(), "nachAsssembleMatrix", "", 7, 0);


        // treat Dirichlet boundary conditions
        incorporateDirichletBC_();
//        std::cout << std::endl;
//        printmatrix(std::cout, this->matrix(), "nachDirichletMatrix", "", 7, 0);
    }

private:

    void incorporateDirichletBC_()
    {
        // the local restrictions of the global fe basis
        auto localView = feBasis().localView();
        auto localIndexSet = feBasis().localIndexSet();

        for (const auto& element : elements(this->gridView_()))
        {
            auto eg = element.geometry();

            localView.bind(element);
            localIndexSet.bind(localView);

            // number of dofs in this element
            auto numLocalDofs = localView.tree().size();

            // loop over the intersections of the element
            for (const auto& is : intersections(this->gridView_(), element))
            {
                // handle only faces on the boundary
                if (!is.boundary())
                    continue;

                // only treat faces with dirichlet boundary conditions
                auto bcTypes = this->problem_().boundaryTypes(element, is);
                if (!bcTypes.hasDirichlet())
                    continue;

                // get face index of this intersection
                int fIdx = is.indexInInside();

                // get reference elements and finite element
                const auto& fe = localView.tree().finiteElement();
                const auto& refElement = ReferenceElements::general(eg.type());

                // handle Dirichlet boundaries
                for (unsigned int i = 0; i < numLocalDofs; i++)
                {
                    const auto dofIdxGlobal = localIndexSet.index(i);

                    const auto& localKey = fe.localCoefficients().localKey(i);
                    auto subEntity = localKey.subEntity();
                    auto codim = localKey.codim();

                    //std::cout << "codim= " << codim << "   " << "dofidxglobal= " << dofIdxGlobal << std::endl;

                    // skip interior dofs
                    if (codim == 0){

                        continue;
                    }
                    // iterate over number of degrees of freedom with the given codim which are on the current face
                    for (int j = 0; j < refElement.size(fIdx, 1, codim); j++)
                    {
                        // If j-th sub entity is the sub entity corresponding to local dof, continue and assign BC
                        if (subEntity == refElement.subEntity(fIdx, 1, j, codim))
                        {
                            // get global coordinate of this degree of freedom
                            auto globalPos = eg.global(refElement.position(subEntity, codim));

                            // value of dirichlet BC
                            const auto diriValues = this->problem_().dirichlet(element, is, globalPos);

//printvector(std::cout, diriValues, "impFemAssemblerDiriValues","");

                            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            {
                                if (bcTypes.isDirichlet(eqIdx))
                                {
                                    // set Dirichlet value on the RHS
                                    auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                                    auto& residual = this->residual()[dofIdxGlobal][eqIdx];
                                    residual = this->model_().curSol()[dofIdxGlobal][eqIdx] - diriValues[pvIdx];
//                           std::cout << "impFemAssemblerPvIdx: " << pvIdx << std::endl;
//                           std::cout << "impFemAssemblerDiriValues[pvIdx]: " << diriValues[pvIdx] << std::endl;
//                           this->residual()[dofIdxGlobal][eqIdx] = this->model_().curSol()[dofIdxGlobal][eqIdx] - diriValues[pvIdx];
//std::cout << residual << std::endl;

//std::cout << std::endl;
//printmatrix(std::cout, this->matrix(), "matrix", "", 7, 0);

                                    // Modify the row of the global matrix
                                    auto cIt = this->matrix()[dofIdxGlobal].begin();
                                    auto cEndIt = this->matrix()[dofIdxGlobal].end();
                                    for (; cIt != cEndIt; ++cIt)
                                    {
                                        (*cIt)[eqIdx] = 0.0;
                                        if (cIt.index() == dofIdxGlobal)
                                            (*cIt)[eqIdx][eqIdx] = 1.0;
                                    }
                                }
                            }

                            // we found the dof
                            break;
                        }
                    }
                }
            }
        }
    }

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        auto n = numDofs();

        // allocate raw matrix
        this->matrix_ = std::make_shared<JacobianMatrix>(n, n, JacobianMatrix::random);

        // A view on the FE basis of a single element
        auto localView = feBasis().localView();
        auto localIndexSet = feBasis().localIndexSet();

        // get occupation pattern of the matrix
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(n, n);

        // Loop over all leaf elements
        for(const auto& element : elements(this->gridView_()))
        {
            // Bind the local FE basis view to the current element
            localView.bind(element);
            localIndexSet.bind(localView);

            // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
            for (std::size_t i = 0; i < localView.tree().size(); i++)
                for (std::size_t j = 0; j < localView.tree().size(); j++)
                    occupationPattern.add(localIndexSet.index(i), localIndexSet.index(j));
        }

        // export occupation pattern to matrix
        occupationPattern.exportIdx(this->matrix());
    }

    std::shared_ptr<FEBasis> feBasisPtr_;
};

} // namespace Dumux

#endif
