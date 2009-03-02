// $Id$

#ifndef DUNE_COUPLEDDIFFUSION_HH
#define DUNE_COUPLEDDIFFUSION_HH

#include "coupledmodel.hh"

namespace Dune
{
/** \todo Please doc me! */

template<class DiffusionModel>
class CoupledDiffusion : public CoupledModel<DiffusionModel, DiffusionModel, CoupledDiffusion<DiffusionModel> > {
public:
    template<int dim>
    struct NodeLayout
    {
        bool contains(GeometryType gt) {
            return gt.dim() == 0;
        }
    };

    // typedefs depending on base class
    typedef CoupledModel<DiffusionModel, DiffusionModel, CoupledDiffusion<DiffusionModel> > BaseType;
    typedef typename BaseType::FirstGrid Grid;
    typedef typename BaseType::FirstMatrixType MatrixType;
    typedef typename BaseType::FirstRowIterator RowIterator;
    typedef typename BaseType::FirstColIterator ColIterator;
    typedef typename BaseType::GlobalMatrixType GlobalMatrixType;
    typedef typename BaseType::GlobalVectorType GlobalVectorType;

    // new typedefs
    enum {dim = Grid::dimension};
    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::template Codim<dim>::Iterator VIterator;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,NodeLayout> VM;
    typedef typename Grid::HostGridType HostGrid;
    typedef typename HostGrid::template Codim<dim>::EntityPointer HostVPointer;

    template <class A12Type, class A21Type>
    void assembleCoupling(A12Type& A12, A21Type& A21)
    {
        // intialize rowsizes with 0
        for (typename A12Type::size_type i = 0; i < A12.N(); i++)
            A12.setrowsize(i, 0);
        for (typename A21Type::size_type i = 0; i < A21.N(); i++)
            A21.setrowsize(i, 0);

        // loop 1 over all nodes of the first grid
        VIterator vendit = (this->firstGrid_).template leafend<dim>();
        for (VIterator firstIt = (this->firstGrid_).template leafbegin<dim>(); firstIt != vendit; ++firstIt)
        {
            // get the node pointer on the host grid
            const HostVPointer& hostVPointer = (this->firstGrid_).template getHostEntity<dim>(*firstIt);

            // check if the node is also part of the second grid
            if (!((this->secondGrid_).template contains<dim>(hostVPointer)))
                continue; // otherwise move on to the next node

            // get the indices of the node with respect to the two local matrices
            int firstId = firstVM.map(*firstIt);
            int secondId = secondVM.map(*firstIt);

            // set nontrivial rowsizes
            A12.setrowsize(firstId, 1);
            A21.setrowsize(secondId, (this->A11)[firstId].size());

            // generate some output
            FieldVector<double,dim> globalCoord = (*firstIt).geometry().corner(0);
        } // end loop 1 over all nodes of the first grid
        A12.endrowsizes();
        A21.endrowsizes();

        // loop 2 over all nodes of the first grid
        for (VIterator firstIt = (this->firstGrid_).template leafbegin<dim>(); firstIt != vendit; ++firstIt)
        {
            // get the node pointer on the host grid
            const HostVPointer& hostVPointer = (this->firstGrid_).template getHostEntity<dim>(*firstIt);

            // check if the node is also part of the second grid
            if (!((this->secondGrid_).template contains<dim>(hostVPointer)))
                continue; // otherwise move on to the next node

            // get the indices of the node with respect to the two local matrices
            int firstId = firstVM.map(*firstIt);
            int secondId = secondVM.map(*firstIt);

            // add the index to A12
            A12.addindex(firstId, secondId);

            // add the indices to A21
            ColIterator endJBlock = (this->A11)[firstId].end();
            for (ColIterator jBlock = (this->A11)[firstId].begin(); jBlock != endJBlock; ++jBlock)
                A21.addindex(secondId, jBlock.index());
        } // end loop 2 over all nodes of the first grid
        A12.endindices();
        A21.endindices();

        // loop 3 over all nodes of the first grid
        for (VIterator firstIt = (this->firstGrid_).template leafbegin<dim>(); firstIt != vendit; ++firstIt)
        {
            // get the node pointer on the host grid
            const HostVPointer& hostVPointer = (this->firstGrid_).template getHostEntity<dim>(*firstIt);

            // check if the node is also part of the second grid
            if (!((this->secondGrid_).template contains<dim>(hostVPointer)))
                continue; // otherwise move on to the next node

            // get the indices of the node with respect to the two local matrices
            int firstId = firstVM.map(*firstIt);
            int secondId = secondVM.map(*firstIt);

            // set the entries
            A12[firstId][secondId] = -1.0;
            A21[secondId] = (this->A11)[firstId];
            (this->A11)[firstId] = 0.0;
            (this->A11)[firstId][firstId] = 1.0;
            (this->secondModel_).rhs()[secondId] += (this->firstModel_).rhs()[firstId];
            (this->firstModel_).rhs()[firstId] = 0;
        } // end loop 3 over all nodes of the first grid

    }

    virtual void solve()
    {
        typedef MatrixAdapter<GlobalMatrixType,GlobalVectorType,GlobalVectorType> Operator;
        Operator op(this->A);
        double red=1E-14;
        SeqILU0<GlobalMatrixType,GlobalVectorType,GlobalVectorType> ilu0(this->A,1.0);
        BiCGSTABSolver<GlobalVectorType> solver(op,ilu0,red,10000,1);
        InverseOperatorResult r;
        solver.apply(this->u, this->f, r);


        if (stationary_)
        {
            this->u *= -1;

            const typename MatrixType::block_type::size_type rowsInBlock1 = MatrixType::block_type::rows;
            const typename MatrixType::block_type::size_type rowsInBlock2 = MatrixType::block_type::rows;
            typename MatrixType::size_type nOfBlockRows1 = (this->firstModel_).matrix().N();

            Iterator dummyIT1((this->firstGrid_).template leafbegin<0>());
            IntersectionIterator dummyIS1(IntersectionIteratorGetter<Grid,LeafTag>::begin(*dummyIT1));
            VIterator endItV1 = (this->firstGrid_).template leafend<dim>();
            for (VIterator it = (this->firstGrid_).template leafbegin<dim>(); it != endItV1; ++it)
            {
                FieldVector<double,dim> globalCoord = (*it).geometry().corner(0);
                FieldVector<BoundaryConditions::Flags, rowsInBlock1> bctype = (this->firstModel_).problem.bctype(globalCoord, *dummyIT1, dummyIS1, globalCoord);
                int firstId = firstVM.map(*it);

                for (typename MatrixType::block_type::size_type i = 0; i < rowsInBlock1; i++)
                    if (bctype[i] == BoundaryConditions::dirichlet)
                    {
                        FieldVector<double, rowsInBlock1> dirichletBC = (this->firstModel_).problem.g(globalCoord, *dummyIT1, dummyIS1, globalCoord);
                        this->u[firstId*rowsInBlock1 + i] += dirichletBC[i];
                    }
            }
            Iterator dummyIT2((this->secondGrid_).template leafbegin<0>());
            IntersectionIterator dummyIS2(IntersectionIteratorGetter<Grid,LeafTag>::begin(*dummyIT2));
            VIterator endItV2 = (this->secondGrid_).template leafend<dim>();
            for (VIterator it = (this->secondGrid_).template leafbegin<dim>(); it != endItV2; ++it)
            {
                FieldVector<double,dim> globalCoord = (*it).geometry().corner(0);
                FieldVector<BoundaryConditions::Flags, rowsInBlock1> bctype = (this->secondModel_).problem.bctype(globalCoord, *dummyIT2, dummyIS2, globalCoord);
                int secondId = secondVM.map(*it);

                for (typename MatrixType::block_type::size_type i = 0; i < rowsInBlock2; i++)
                    if (bctype[i] == BoundaryConditions::dirichlet)
                    {
                        FieldVector<double, rowsInBlock2> dirichletBC = (this->secondModel_).problem.g(globalCoord, *dummyIT2, dummyIS2, globalCoord);
                        this->u[rowsInBlock1*nOfBlockRows1 + secondId*rowsInBlock2 + i] += dirichletBC[i];
                    }
            }
        }
        // transfer to local solution vectors
        const typename MatrixType::block_type::size_type colsInBlock1 = MatrixType::block_type::cols;
        const typename MatrixType::block_type::size_type colsInBlock2 = MatrixType::block_type::cols;
        for (int i = 0; i < (this->firstModel_).sol().size(); i++)
            for (typename MatrixType::block_type::size_type k = 0; k < colsInBlock1; k++)
                (this->firstModel_).sol()[i][k] = this->u[i*colsInBlock1 + k];
        for (int i = 0; i < (this->secondModel_).sol().size(); i++)
            for (typename MatrixType::block_type::size_type k = 0; k < colsInBlock2; k++)
                (this->secondModel_).sol()[i][k] = this->u[colsInBlock1*(this->firstModel_).sol().size() + i*colsInBlock2 + k];
    }


    CoupledDiffusion(const Grid& firstGrid, DiffusionModel& firstModel,
                     const Grid& secondGrid, DiffusionModel& secondModel,
                     bool assembleGlobalSystem, bool stationary)
        : BaseType(firstGrid, firstModel, secondGrid, secondModel, assembleGlobalSystem),
          firstVM(firstGrid, firstGrid.leafIndexSet()), secondVM(secondGrid, secondGrid.leafIndexSet()),
          stationary_(stationary)
    {}

protected:
    VM firstVM;
    VM secondVM;
    bool stationary_;

};

}

#endif
