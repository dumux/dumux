// $Id$

#ifndef DUNE_LOCALMIXED_HH
#define DUNE_LOCALMIXED_HH


#include<dune/disc/operators/localstiffness.hh>
#include"dumux/stokes/stokesproblem.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar, int m>
class LocalMixed : public LinearLocalStiffness<typename Grid::LeafGridView, Scalar, m>
{
    //    typedef LocalStiffness< LocalMixed<Grid,Scalar,m>, Grid, Scalar, m> Stiffness;
    typedef LinearLocalStiffness<typename Grid::LeafGridView, Scalar, m> Stiffness;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    typedef StokesProblem<Grid, Scalar> Problem;
    enum {dim=Grid::dimension};

public:
    LocalMixed(Problem& problem)
        : Stiffness(), problem_(problem)
    {}

    // ASSUME a 2D uniform rectangular axiparallel grid
    void assemble (const Element& element, int k = 1)
    {
        // get the number of element faces:
        int nFaces = element.template count<1>();

        this->setcurrentsize(nFaces+1);

        // initialize everything with 0:
        for (int i = 0; i < nFaces+1; i++) {
            this->bctype[i].assign(BoundaryConditions::neumann);
            this->b[i] = 0;
            for (int j = 0; j < nFaces+1; j++)
                this->A[i][j] = 0;
        }

        // obtain the inverse of the Jacobian of the element transformation:
        const Geometry& geometry = element.geometry();
        GeometryType geomType = geometry.type();
        const FieldVector<Scalar,dim>& local = ReferenceElements<Scalar,dim>::general(geomType).position(0,0);
        const FieldMatrix<Scalar,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(local);

        // extract the discretization parameters:
        Scalar deltaX = 1.0/jacobianInv[0][0];
        Scalar deltaY = 1.0/jacobianInv[1][1];
        Scalar volume = deltaX*deltaY;

        // mass conservation
        // control volume: element
        // div v = 0 -> \sum_e |e| v.n = 0
        this->A[4][0] = deltaY;
        this->A[4][1] = -deltaY;
        this->A[4][2] = deltaX;
        this->A[4][3] = -deltaX;

        // momentum conservation
        // - div T + grad p
        // x-comp control volume: element shifted by deltaX/2 in x-direction
        // y-comp control volume: element shifted by deltaY/2 in y-direction
        // -> i-comp: \sum_e |e|( -T_i.n + p.n_i )
        Scalar mu = 1.0;
        this->A[0][0] = this->A[1][1] = 2.0*mu*deltaY/deltaX;
        this->A[2][2] = this->A[3][3] = 2.0*mu*deltaX/deltaY;
        this->A[0][1] = this->A[1][0] = -2.0*mu*deltaY/deltaX;
        this->A[2][3] = this->A[3][2] = -2.0*mu*deltaX/deltaY;
        this->A[0][4] = deltaY;
        this->A[1][4] = -deltaY;
        this->A[2][4] = deltaX;
        this->A[3][4] = -deltaX;

        // insert the source term and the boundary conditions:
        IntersectionIterator endIsIt = element.ileafend();
        for (IntersectionIterator isIt = element.ileafbegin(); isIt != endIsIt; ++isIt)
        {
            // the local face number:
            int indexInInside = isIt->indexInInside();

            // geometry data of the face:
            GeometryType geomType = isIt->geometryInInside().type();
            const FieldVector<Scalar,dim-1>& localDimM1 = ReferenceElements<Scalar,dim-1>::general(geomType).position(0,0);
            const FieldVector<Scalar,dim>& local = ReferenceElements<Scalar,dim>::general(geomType).position(indexInInside, 1);
            FieldVector<Scalar,dim> global = isIt->geometry().global(localDimM1);

            // get the source term from the problem:
            FieldVector<Scalar,dim+1> source = problem_.q(global, element, local);

            // set the right hand side:
            this->b[indexInInside] = 0.5*volume*source[indexInInside/2];

            // after this, only boundary elements should be considered:
            if (isIt->neighbor())
                continue;

            // get the boundary condition type:
            BoundaryConditions::Flags bctype = problem_.bctype(global, element, isIt, local);

            if (bctype == BoundaryConditions::dirichlet)
            {
                // get the boundary value from the problem:
                FieldVector<Scalar,dim> velocity = problem_.dirichlet(global, element, isIt, local);

                // set the right hand side and the boundary condition type:
                this->b[indexInInside] = velocity[indexInInside/2];
                this->bctype[indexInInside].assign(BoundaryConditions::dirichlet);
            }
            else
            {
                DUNE_THROW(NotImplemented, "BoundaryConditions != Dirichlet for LocalMixed");
            }
        }

        return;
    }

    void assembleBoundaryCondition (const Element& element, int k = 1)
    {
    }

private:
    Problem& problem_;
};

}
#endif
