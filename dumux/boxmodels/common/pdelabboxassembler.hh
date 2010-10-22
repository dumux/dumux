// $Id$
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
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
/*!
 * \file
 *
 * \brief An assembler for the Jacobian matrix based on PDELab.
 */
#ifndef DUMUX_PDELAB_BOX_ASSEMBLER_HH
#define DUMUX_PDELAB_BOX_ASSEMBLER_HH

#include "pdelabboxlocaloperator.hh"

namespace Dumux {

namespace PDELab {

template<class TypeTag>
class BoxAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Constraints)) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ScalarGridFunctionSpace)) ScalarGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalOperator)) LocalOperator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridOperatorSpace)) GridOperatorSpace;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    enum{dim = GridView::dimension};
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef SolutionVector Vector;
    typedef JacobianMatrix Matrix;
    typedef Matrix RepresentationType;

    enum {
        enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)),
        enableJacobianRecycling = GET_PROP_VALUE(TypeTag, PTAG(EnableJacobianRecycling)),

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))
    };

    // copying the jacobian assembler is not a good idea
    BoxAssembler(const BoxAssembler &);

public:
    enum EntityColor {
        Red, //!< entity needs to be reassembled because it's above the tolerance
        Yellow, //!< entity needs to be reassembled because a neighboring element is red
        Green //!< entity does not need to be reassembled
    };

    BoxAssembler()
    {
        problemPtr_ = 0;
        fem_ = 0;
        cn_ = 0;
        scalarGridFunctionSpace_ = 0;
        gridFunctionSpace_ = 0;
        constraintsTrafo_ = 0;
        localOperator_ = 0;
        gridOperatorSpace_ = 0;
        matrix_ = 0;
    }

    ~BoxAssembler()
    {
        delete matrix_;
        delete gridOperatorSpace_;
        delete localOperator_;
        delete constraintsTrafo_;
        delete gridFunctionSpace_;
        delete scalarGridFunctionSpace_;
        delete cn_;
        delete fem_;
    }

    void init(Problem& problem)
    {
        inJacobianAssemble_ = false;
        problemPtr_ = &problem;
        fem_ = new FEM();
        //cn_ = new Constraints(*problemPtr_);
        cn_ = new Constraints();
        scalarGridFunctionSpace_ = new ScalarGridFunctionSpace(problemPtr_->gridView(), *fem_, *cn_);
        gridFunctionSpace_ = new GridFunctionSpace(*scalarGridFunctionSpace_);

        //cn_->compute_ghosts(*gridFunctionSpace_);

        //typedef BoundaryIndexHelper<TypeTag> BoundaryFunction;
        //BoundaryFunction *bTypes = new BoundaryFunction();
        constraintsTrafo_ = new ConstraintsTrafo();
        //Dune::PDELab::constraints(*bTypes, *gridFunctionSpace_, *constraintsTrafo_, false);

        // initialize the grid operator spaces
        localOperator_ = new LocalOperator(problemPtr_->model());
        gridOperatorSpace_ =
            new GridOperatorSpace(*gridFunctionSpace_, *constraintsTrafo_,
                                  *gridFunctionSpace_, *constraintsTrafo_, *localOperator_);

        // initialize the jacobian matrix and the right hand side
        // vector
        matrix_ = new Matrix(*gridOperatorSpace_);
        *matrix_ = 0;
        reuseMatrix_ = false;

        // calculate the ghost vertices
        VertexIterator vIt = gridView_().template begin<dim>();
        VertexIterator vEndIt = gridView_().template end<dim>();
        for (; vIt != vEndIt; ++vIt) {
            if (vIt->partitionType() == Dune::GhostEntity) {
                int vIdx = vertexMapper_().map(*vIt);
                ghostIndices_.push_back(vIdx);
            }
        };

        int numVerts = gridView_().size(dim);
        int numElems = gridView_().size(0);
        residual_.resize(numVerts);
      
        // initialize data needed for partial reassembly
        if (enablePartialReassemble) {
            evalPoint_.resize(numVerts);
            vertexColor_.resize(numVerts);
            elementColor_.resize(numElems);
        }
        reassembleAll();
    }

    SolutionVector &evalPoint()
    { return evalPoint_; }

    const SolutionVector &evalPoint() const
    { return evalPoint_; }

    bool inJacobianAssemble() const
    { return inJacobianAssemble_; }
    
    void assemble()
    {
        // assemble the global jacobian matrix
        if (!reuseMatrix_) {
            if (enablePartialReassemble) {
                // move the evaluation points of red vertices
                for (int i = 0; i < vertexColor_.size(); ++i) {
                    if (vertexColor_[i] == Red)
                        evalPoint_[i] = model_().curSol()[i];
                }
            }

            reassembleTolerance_ = nextReassembleTolerance_;

            // we actually need to reassemle!
            resetMatrix_();
            inJacobianAssemble_ = true;
            gridOperatorSpace_->jacobian(model_().curSol(), *matrix_);
            inJacobianAssemble_ = false;
        }
        reuseMatrix_ = false;

        // calculate the global residual
        residual_ = 0;
        gridOperatorSpace_->residual(model_().curSol(), residual_);

        typedef typename Matrix::block_type BlockType;
        // set the entries for the ghost nodes
        BlockType Id(0.0);
        for (int i=0; i < BlockType::rows; ++i)
            Id[i][i] = 1.0;

        for (int i=0; i < ghostIndices_.size(); ++i) {
            int globI = ghostIndices_[i];

            (*matrix_)[globI] = 0;
            (*matrix_)[globI][globI] = Id;
            residual_[globI] = 0;
            model_().curSol()[globI] = 0;
        }
    }

    void setMatrixReuseable(bool yesno = true)
    {
        if (enableJacobianRecycling)
            reuseMatrix_ = yesno;
    }

    
    void reassembleAll()
    {
        nextReassembleTolerance_ = 0.0;
        if (enablePartialReassemble) {
            std::fill(vertexColor_.begin(),
                      vertexColor_.end(),
                      Red);
            std::fill(elementColor_.begin(),
                      elementColor_.end(), 
                      Red);
        }
    }
   
    Scalar reassembleTolerance() const
    { return reassembleTolerance_; }

    void markVertexRed(int globalVertIdx, bool yesno)
    {
        if (enablePartialReassemble)
            vertexColor_[globalVertIdx] = yesno?Red:Green;
    }
    
    void computeColors(Scalar relTol)
    { 
        if (!enablePartialReassemble)
            return;

        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();

        nextReassembleTolerance_ = 0;

        int redVert = 0;
        // mark the red vertices
        for (int i = 0; i < vertexColor_.size(); ++i) {
            const PrimaryVariables &evalPv = evalPoint_[i];
            const PrimaryVariables &solPv = model_().curSol()[i];
            
            Scalar vertErr = model_().relativeErrorVertex(i,
                                                          evalPv,
                                                          solPv);

            // mark vertex as red or green
            vertexColor_[i] = (vertErr > relTol)?Red:Green;
            if (vertErr > relTol)
                ++redVert;
            else
                nextReassembleTolerance_ = 
                    std::max(nextReassembleTolerance_, vertErr);
        };

        // Mark all red elements
        for (; elemIt != elemEndIt; ++elemIt) {
            bool needReassemble = false;
            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                if (vertexColor_[globalI] == Red) {
                    needReassemble = true;
                    break;
                }
            };
            
            int globalElemIdx = elementMapper_().map(*elemIt);
            elementColor_[globalElemIdx] = needReassemble?Red:Green;
        }
        
        // Mark all yellow vertices
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] == Green)
                continue; // green elements do not tint vertices
                          // yellow!
            
            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                // if a vertex is already red, don't recolor it to
                // yellow!
                if (vertexColor_[globalI] != Red)
                    vertexColor_[globalI] = Yellow;
            };
        }

        // Mark all yellow elements
        int numGreen = 0;
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] == Red) {
                continue; // element is red already!
            }

            bool isYellow = false;
            int numVerts = elemIt->template count<dim>();
            for (int i=0; i < numVerts; ++i) {
                int globalI = vertexMapper_().map(*elemIt, i, dim);
                if (vertexColor_[globalI] == Yellow) {
                    isYellow = true;
                    break;
                }
            };
            
            if (isYellow) {
                elementColor_[elemIdx] = Yellow;
            }
            else // elementColor_[elemIdx] == Green
                ++ numGreen;
        }
        
        int numTot = gridView_().size(0);
        numTot = gridView_().comm().sum(numTot);
        numGreen = gridView_().comm().sum(numGreen);
        elemsReasm = numTot - numGreen;
        problem_().newtonController().endIterMsg()
            << ", reassemble " 
            << numTot - numGreen << "/" << numTot
            << " (" << 100*Scalar(numTot - numGreen)/numTot << "%) elems";
    };
    
    int elemsReasm;
    int vertexColor(const Element &element, int vertIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!
        
        int globalIdx = vertexMapper_().map(element, vertIdx, dim);
        return vertexColor_[globalIdx];
    }

    int vertexColor(int globalVertIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!
        return vertexColor_[globalVertIdx];
    }

    int elementColor(const Element &element) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!
        
        int globalIdx = elementMapper_().map(element);
        return elementColor_[globalIdx];
    }

    int elementColor(int globalElementIdx) const
    {
        if (!enablePartialReassemble)
            return Red; // reassemble unconditionally!
        return elementColor_[globalElementIdx];
    }
   
    const GridFunctionSpace& gridFunctionSpace() const
    {
        return *gridFunctionSpace_;
    }

    const ConstraintsTrafo& constraintsTrafo() const
    {
        return *constraintsTrafo_;
    }

    //! return const reference to matrix
    const Matrix& matrix() const
    { return *matrix_; }

    //! return const reference to residual
    const SolutionVector& residual() const
    { return residual_; }


private:
    void resetMatrix_()
    {
        if (!enablePartialReassemble) {
            (*matrix_) = 0;
            return;
        }
      
        // reset all entries corrosponding to a red vertex
        for (int rowIdx = 0; rowIdx < matrix_->N(); ++rowIdx) {
            if (vertexColor_[rowIdx] == Green)
                continue; // the equations for this control volume are
                          // already below the treshold
            // reset row to 0
            typedef typename JacobianMatrix::ColIterator ColIterator;
            ColIterator colIt = (*matrix_)[rowIdx].begin();
            const ColIterator &colEndIt = (*matrix_)[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                (*colIt) = 0.0;
            }
        };

        //printSparseMatrix(std::cout, *matrix_, "J", "row");
    }
    
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
    Constraints *cn_;
    FEM *fem_;
    ScalarGridFunctionSpace *scalarGridFunctionSpace_;
    GridFunctionSpace *gridFunctionSpace_;
    ConstraintsTrafo *constraintsTrafo_;
    LocalOperator *localOperator_;
    GridOperatorSpace *gridOperatorSpace_;

    Matrix *matrix_;
    bool reuseMatrix_;
    bool inJacobianAssemble_;
    SolutionVector evalPoint_;
    std::vector<EntityColor> vertexColor_;
    std::vector<EntityColor> elementColor_;
    std::vector<int> ghostIndices_;

    Scalar nextReassembleTolerance_;
    Scalar reassembleTolerance_;
    SolutionVector residual_;
};

} // namespace PDELab
} // namespace Dumux

#endif
