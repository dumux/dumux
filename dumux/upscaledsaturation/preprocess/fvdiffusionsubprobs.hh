// $Id$

#ifndef DUNE_FVDIFFUSIONSUBPROBS_HH
#define DUNE_FVDIFFUSIONSUBPROBS_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/upscaledsaturation/preprocess/diffusionsubprobs.hh"
#include "dumux/pardiso/pardiso.hh"
//#include "dumux/diffusion/problems/uniformproblem.hh"
//#include "dumux/transport/problems/simpleproblem.hh"

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Bernd Flemisch, Jochen Fritz; last changed by Markus Wolff
 */

namespace Dune
{
//! \ingroup diffusion
//! Finite Volume Diffusion Model
/*! Provides a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and
 * \f$\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation, \f$q\f$ the source term.
 Template parameters are:

 - Grid         a DUNE grid type
 - Scalar        type used for return values
*/
template<class Grid, class Scalar, class VC> class FVDiffSubProbs: public DiffusionSubProbs<
    Grid, Scalar, VC>
{
    template<int dim> struct ElementLayout
    {
        bool contains(GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
        {
            dim = Grid::dimension
        };
    enum
        {
            dimWorld = Grid::dimensionworld
        };

    typedef    typename Grid::HostGridType HostGrid;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> ElementMapper;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;
    typedef BCRSMatrix<MB> MatrixType;
    typedef BlockVector<FieldVector<Scalar, 1> > Vector;

public:

    void assemble(const Scalar t);

    void solve();

    void pressure(const Scalar t=0)
    {
        assemble(t);
        solve();
        return;
    }

    void calcTotalVelocity(const Scalar t) const;

    void initializeMatrix();

    FVDiffSubProbs(Grid& grid, FractionalFlowProblemSubProbs<Grid, Scalar, VC>& problem) :
        DiffusionSubProbs<Grid, Scalar, VC>(grid, problem),
        elementMapper(grid, grid.levelIndexSet(this->level())), A(grid.size(
                                                                            this->level(), 0), grid.size(this->level(), 0), (2*dim+1)
                                                                  *grid.size(this->level(), 0), BCRSMatrix<MB>::random),
        f(grid.size(this->level(), 0)), solverName_("BiCGSTAB"),
        preconditionerName_("SeqILU0")
    {
        initializeMatrix();
    }

    ElementMapper elementMapper;

private:

    MatrixType A;
    BlockVector< FieldVector<Scalar,1> > f;
    std::string solverName_;
    std::string preconditionerName_;
};

template<class Grid, class Scalar, class VC> void FVDiffSubProbs<Grid, Scalar, VC>::initializeMatrix()
{

    const GridView& gridView = this->grid.levelView(this->level());

    // determine matrix row sizes
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator
            isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                 isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
            if (isIt->neighbor())
                rowSize++;
        A.setrowsize(globalIdxI, rowSize);
    }
    A.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        // add diagonal index
        A.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator
            isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                 isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = elementMapper.map(*outside);

                // add off diagonal index
                A.addindex(globalIdxI, globalIdxJ);
            }
    }
    A.endindices();

    return;
}

template<class Grid, class Scalar, class VC> void FVDiffSubProbs<Grid, Scalar, VC>::assemble(const Scalar t=0)
{
    // initialization: set matrix A to zero
    A = 0;

    const GridView& gridView = this->grid.levelView(this->level());

    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        const HostElementPointer& eItHost = (this->grid).template getHostEntity<0>(*eIt);

        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)*ReferenceElements<Scalar,dim>::general(gt).volume();

        // set right side to zero
        f[globalIdxI] = volume*this->diffProblem.sourcePress(globalPos, *eIt, localPos);

        // get absolute permeability
        FieldMatrix Ki(this->diffProblem.soil.K(globalPos, *eItHost, localPos));

        //compute total mobility
        Scalar lambdaI;
        Scalar satI = this->diffProblem.variables.saturation[globalIdxI];
        lambdaI = this->diffProblem.materialLaw.mobTotal(satI,globalPos, *eItHost, localPos);

        IntersectionIterator
            isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                 isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
        {

            // get geometry type of face
            GeometryType faceGT = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const FieldVector<Scalar,dim-1>&
                faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition& localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(isIt->numberInSelf(),1);

            // get normal vector
            FieldVector<Scalar,dimWorld> unitOuterNormal
                = isIt->unitOuterNormal(faceLocal);

            // get normal vector scaled with volume
            FieldVector<Scalar,dimWorld> integrationOuterNormal= isIt->integrationOuterNormal(faceLocal);
            integrationOuterNormal*= ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

            // get face volume
            Scalar faceVol = 1;
            switch (Grid::dimension)
            {
            case 1: break;
            default: faceVol = isIt->intersectionGlobal().volume();
                break;
            }

            // compute directed permeability vector Ki.n
            FieldVector<Scalar,dim> Kni(0);
            Ki.umv(unitOuterNormal, Kni);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = elementMapper.map(*neighborPointer);

                // compute factor in neighbor
                GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                const HostElementPointer& neighborPointerHost = (this->grid).template getHostEntity<0>(*neighborPointer);

                // distance vector between barycenters
                FieldVector<Scalar,dimWorld>
                    distVec = globalPos - globalPosNeighbor;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // get absolute permeability
                FieldMatrix Kj(this->diffProblem.soil.K(globalPosNeighbor, *neighborPointerHost, localPosNeighbor));

                // compute vectorized permeabilities
                FieldVector<Scalar,dim> Knj(0);
                Kj.umv(unitOuterNormal, Knj);
                Scalar K_n_i = Kni * unitOuterNormal;
                Scalar K_n_j = Knj * unitOuterNormal;
                Scalar Kn = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                // compute permeability tangential to intersection and take arithmetic mean
                FieldVector<Scalar,dim> uON = unitOuterNormal;
                FieldVector<Scalar,dim> K_t_i = Kni - (uON *= K_n_i);
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> K_t_j = Knj - (uON *= K_n_j);
                FieldVector<Scalar,dim> Kt = (K_t_i += K_t_j);
                Kt *= 0.5;
                // Build vectorized averaged permeability
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> K = (Kt += (uON *=Kn));

                //compute total mobility
                Scalar lambdaJ;
                Scalar satJ = this->diffProblem.variables.saturation[globalIdxJ];

                lambdaJ = this->diffProblem.materialLaw.mobTotal(satJ,globalPosNeighbor, *neighborPointerHost, localPosNeighbor);

                // compute averaged total mobility
                // CAREFUL: Harmonic weightig can generate zero matrix entries,
                // use arithmetic weighting instead:

                Scalar lambda = 0.5*(lambdaI + lambdaJ);

                // update diagonal entry
                Scalar entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
                A[globalIdxI][globalIdxI] += entry;

                // set off-diagonal entry
                A[globalIdxI][globalIdxJ] = -entry;
            }
            // boundary face

            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->intersectionGlobal().global(faceLocal);

                // compute total mobility
                Scalar lambda = lambdaI;

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);
                if (bctype == BoundaryConditions::dirichlet)
                {
                    FieldVector<Scalar,dimWorld> distVec(globalPos - globalPosFace);
                    Scalar dist = distVec.two_norm();
                    A[globalIdxI][globalIdxI] -= lambda*faceVol*(Kni*distVec)/(dist*dist);
                    Scalar g = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    f[globalIdxI] -= lambda*faceVol*g*(Kni*distVec)/(dist*dist);
                }
                else
                {
                    Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    f[globalIdxI] -= faceVol*J;
                }
            }
        }
        // end all intersections
    } // end grid traversal
    return;
}

template<class Grid, class Scalar, class VC> void FVDiffSubProbs<Grid, Scalar, VC>::solve()
{
    MatrixAdapter<MatrixType,Vector,Vector> op(A);
    InverseOperatorResult r;

    if (preconditionerName_ == "SeqILU0")
    {
        SeqILU0<MatrixType,Vector,Vector> preconditioner(A, 1.0);
        if (solverName_ == "CG")
        {
            CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
            solver.apply((this->diffProblem.variables.pressure), f, r);
        }
        else if (solverName_ == "BiCGSTAB")
        {
            BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply((this->diffProblem.variables.pressure), f, r);
        }
        else
            DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                       << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
        SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
        if (solverName_ == "Loop")
        {
            LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(this->diffProblem.variables.pressure, f, r);
        }
        else
            DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                       << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner "
                   << preconditionerName_ << ".");
    //    printmatrix(std::cout, A, "global stiffness matrix", "row", 11, 3);
    //    printvector(std::cout, f, "right hand side", "row", 200, 1, 3);
    //    printvector(std::cout, (this->diffProblem.variables.pressure), "pressure", "row", 200, 1, 3);
    return;
}

template<class Grid, class Scalar, class VC> void FVDiffSubProbs<Grid, Scalar, VC>::calcTotalVelocity(const Scalar t=0) const
{
    const GridView& gridView = this->grid.levelView(this->level());

    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt
             != eItEnd; ++eIt)
    {
        // cell geometry type
        GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition&
            localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // cell center in global coordinates
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        const HostElementPointer& eItHost = (this->grid).template getHostEntity<0>(*eIt);

        // cell index
        int globalIdxI = this->elementMapper.map(*eIt);

        // get pressure and permeability in element
        Scalar pressI = this->diffProblem.variables.pressure[globalIdxI];

        // get absolute permeability
        FieldMatrix kI(this->diffProblem.soil.K(globalPos, *eItHost, localPos));

        //compute total mobility
        Scalar lambdaI;
        Scalar satI = this->diffProblem.variables.saturation[globalIdxI];
        lambdaI = this->diffProblem.materialLaw.mobTotal(satI,globalPos, *eItHost, localPos);

        Scalar faceVol[2*dim];

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt!=isItEnd; ++isIt)
        {
            // get geometry type of face
            GeometryType faceGT = isIt->intersectionSelfLocal().type();

            //Geometry dg = isIt->intersectionSelfLocal();
            // local number of face
            int numberInSelf = isIt->numberInSelf();

            switch (dim)
            {
            case 1:
                faceVol[numberInSelf] = 1;
            default:
                faceVol[numberInSelf] = isIt->intersectionGlobal().volume();
            }

            // center in face's reference element
            const FieldVector<Scalar,dim-1>&
                faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition&
                localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(numberInSelf,1);

            // get normal vector
            FieldVector<Scalar,dimWorld> unitOuterNormal
                = isIt->unitOuterNormal(faceLocal);

            // center of face in globalPos coordinates
            const GlobalPosition&
                globalPosFace = isIt->intersectionGlobal().global(faceLocal);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->elementMapper.map(*neighborPointer);

                // get neighbor pressure and permeability
                Scalar pressJ = this->diffProblem.variables.pressure[globalIdxJ];

                // compute factor in neighbor
                GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition&
                    localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                // neighbor cell center in globalPos coordinates
                const GlobalPosition&
                    globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                const HostElementPointer& neighborPointerHost = (this->grid).template getHostEntity<0>(*neighborPointer);

                // distance vector between barycenters
                FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosNeighbor;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // get absolute permeability
                FieldMatrix kJ(this->diffProblem.soil.K(globalPosNeighbor, *neighborPointerHost, localPosNeighbor));

                // compute vectorized permeabilities
                FieldVector<Scalar,dim> kNI(0);
                FieldVector<Scalar,dim> kNJ(0);
                kI.umv(unitOuterNormal, kNI);
                kJ.umv(unitOuterNormal, kNJ);
                // compute permeability normal to intersection and take harmonic mean
                Scalar k_N_I = kNI * unitOuterNormal;
                Scalar k_N_J = kNJ * unitOuterNormal;
                Scalar kN = 2 * k_N_I * k_N_J / (k_N_I + k_N_J);
                // compute permeability tangential to intersection and take arithmetic mean
                FieldVector<Scalar,dim> uON = unitOuterNormal;
                FieldVector<Scalar,dim> k_T_I = kNI - (uON *= k_N_I);
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> k_T_J = kNJ - (uON *= k_N_J);
                FieldVector<Scalar,dim> kT = (k_T_I += k_T_J);
                kT *= 0.5;
                // Build vectorized averaged permeability
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> permeability = (kT += (uON *=kN));

                //compute total mobility
                Scalar lambdaJ;
                Scalar satJ = this->diffProblem.variables.saturation[globalIdxJ];
                lambdaJ = this->diffProblem.materialLaw.mobTotal(satJ,globalPosNeighbor, *neighborPointerHost, localPosNeighbor);

                // compute averaged total mobility
                // CAREFUL: Harmonic weightig can generate zero matrix entries,
                // use arithmetic weighting instead:
                Scalar lambda = 1;

                lambda = 0.5*(lambdaI + lambdaJ);

                FieldVector<Scalar,dimWorld> vTotal(permeability);
                vTotal *= lambda*(pressI - pressJ)/dist;

                this->diffProblem.variables.velocity[globalIdxI][numberInSelf] = vTotal;
            }
            // boundary face

            else
            {
                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);
                if (bctype == BoundaryConditions::dirichlet)
                {
                    // distance vector between barycenters
                    FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;

                    Scalar dist = distVec.two_norm();
                    distVec /= dist;

                    // compute directed permeability vector kI.n
                    FieldVector<Scalar,dim> kNI(0);
                    kI.umv(distVec, kNI);

                    // compute averaged total mobility
                    Scalar lambda = 1.;
                    lambda = lambdaI;

                    Scalar g = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);

                    FieldVector<Scalar,dim> vTotal(kNI);
                    vTotal *= lambda*(g-pressI)/dist;

                    this->diffProblem.variables.velocity[globalIdxI][numberInSelf] = vTotal;
                }
                else
                {
                    Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    FieldVector<Scalar,dimWorld> unitOuterNormal
                        = isIt->unitOuterNormal(faceLocal);
                    this->diffProblem.variables.velocity[globalIdxI][numberInSelf] = unitOuterNormal;
                    this->diffProblem.variables.velocity[globalIdxI][numberInSelf] *= J;
                }

            }
        }
        //        std::cout<<"velocity = "<< this->diffProblem.variables.velocity <<std::endl;
        if (dim == 2)
        {
            Scalar sum = (fabs(this->diffProblem.variables.velocity[globalIdxI][0][0]*faceVol[0])
                          + fabs(this->diffProblem.variables.velocity[globalIdxI][1][0]*faceVol[1])
                          + fabs(this->diffProblem.variables.velocity[globalIdxI][2][1]*faceVol[2])
                          + fabs(this->diffProblem.variables.velocity[globalIdxI][3][1]*faceVol[3]));
            Scalar diff = fabs(this->diffProblem.variables.velocity[globalIdxI][0][0]*faceVol[0]
                               - this->diffProblem.variables.velocity[globalIdxI][1][0]*faceVol[1]
                               + this->diffProblem.variables.velocity[globalIdxI][2][1]*faceVol[2]
                               - this->diffProblem.variables.velocity[globalIdxI][3][1]*faceVol[3])/sum;
            if (diff> 1e-6&& sum> 1e-9)
            {
                std::cout << "NOT conservative!!! diff = "<< diff
                          << ", globalIdxI = "<< globalIdxI << std::endl;
                std::cout << this->diffProblem.variables.velocity[globalIdxI][0][0]*faceVol[0]<< ", "
                          << this->diffProblem.variables.velocity[globalIdxI][1][0]*faceVol[1]<< ", "
                          << this->diffProblem.variables.velocity[globalIdxI][2][1]*faceVol[2]<< ", "
                          << this->diffProblem.variables.velocity[globalIdxI][3][1]*faceVol[3]<< std::endl;
            }
        }
    } // end grid traversal
    //    std::cout<<this->diffProblem.variables.velocity<<std::endl;
    return;
}

}
#endif
