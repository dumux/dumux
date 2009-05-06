// $Id$

#ifndef DUNE_DIFFUSIONVELOCITYPROBLEM_HH
#define DUNE_DIFFUSIONVELOCITYPROBLEM_HH

#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
/** \todo Please doc me! */

template<class Grid, class Scalar, class VC, class Problem = DiffusionProblem<Grid, Scalar, VC> > class FVDiffusionVelocity: public FVDiffusion<
    Grid, Scalar, VC, Problem>
{

    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum
        {   dim = Grid::dimension};
    enum
        {   dimWorld = Grid::dimensionworld};

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    FVDiffusionVelocity(Grid& grid, Problem& problem)
        : FVDiffusion<Grid,Scalar,VC, Problem>(grid, problem)
    {}

    FVDiffusionVelocity(Grid& grid, Problem& problem, std::string solverName,
                        std::string preconditionerName)
        : FVDiffusion<Grid,Scalar,VC, Problem>(grid, problem,solverName,preconditionerName)
    {}

    void calcTotalVelocity(const Scalar t=0) const
    {
        const GridView& gridView = this->grid.levelView(this->level());

        // find out whether gravity effects are relevant
        bool hasGravity = false;
        const FieldVector<Scalar,dim>& gravity(this->diffProblem.gravity());
        for (int k = 0; k < dim; k++)
            if (gravity[k] != 0)
                hasGravity = true;

        ElementIterator eItEnd = gridView.template end<0>();
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {

            // cell geometry type
            GeometryType gt = eIt->geometry().type();

            // cell center in reference element
            const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

            // cell center in global coordinates
            const GlobalPosition& globalPos = eIt->geometry().global(localPos);

            // cell index
            int globalIdxI = this->diffProblem.variables.diffMapper.map(*eIt);

            // get pressure and permeability in element
            Scalar pressI = this->diffProblem.variables.pressure[globalIdxI];

            // get absolute permeability
            FieldMatrix Ki(this->diffProblem.soil().K(globalPos, *eIt, localPos));

            //compute total mobility
            Scalar lambdaI, fractionalWI;
            Scalar satI = this->diffProblem.variables.saturation[globalIdxI];
            lambdaI = this->diffProblem.materialLaw().mobTotal(satI,globalPos, *eIt, localPos);
            if (hasGravity)
                fractionalWI = this->diffProblem.materialLaw().fractionalW(satI,globalPos, *eIt, localPos);

            Scalar faceVol[2*dim];

            // run through all intersections with neighbors and boundary
            IntersectionIterator isItEnd = gridView.template iend(*eIt);
            for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt!=isItEnd; ++isIt)
            {
                // get geometry type of face
                GeometryType faceGT = isIt->geometryInInside().type();

                //Geometry dg = isIt->geometryInInside();
                // local number of facet
                int numberInSelf = isIt->indexInInside();

                faceVol[numberInSelf] =  isIt->geometry().volume();

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
                    globalPosFace = isIt->geometry().global(faceLocal);

                // handle interior face
                if (isIt->neighbor())
                {
                    // access neighbor
                    ElementPointer neighborPointer = isIt->outside();
                    int globalIdxJ = this->diffProblem.variables.diffMapper.map(*neighborPointer);

                    // get neighbor pressure and permeability
                    Scalar pressJ = this->diffProblem.variables.pressure[globalIdxJ];

                    // compute factor in neighbor
                    GeometryType neighborGT = neighborPointer->geometry().type();
                    const LocalPosition&
                        localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                    // neighbor cell center in globalPos coordinates
                    const GlobalPosition&
                        globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                    // distance vector between barycenters
                    FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosNeighbor;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // get absolute permeability
                    FieldMatrix Kj(this->diffProblem.soil().K(globalPosNeighbor, *neighborPointer, localPosNeighbor));

                    // compute vectorized permeabilities
                    FieldVector<Scalar,dim> Kni(0);
                    FieldVector<Scalar,dim> Knj(0);
                    Ki.umv(unitOuterNormal, Kni);
                    Kj.umv(unitOuterNormal, Knj);
                    // compute permeability normal to intersection and take harmonic mean
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
                    Scalar lambdaJ, fractionalWJ;
                    Scalar satJ = this->diffProblem.variables.saturation[globalIdxJ];
                    lambdaJ = this->diffProblem.materialLaw().mobTotal(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);
                    if (hasGravity)
                        fractionalWJ = this->diffProblem.materialLaw().fractionalW(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);

                    // compute averaged total mobility
                    // CAREFUL: Harmonic weighting can generate zero matrix entries,
                    // use arithmetic weighting instead:
                    Scalar lambda = 1;
                    Scalar fractionalW;
                    lambda = 0.5*(lambdaI + lambdaJ);
                    if (hasGravity)
                        fractionalW = 0.5*(fractionalWI + fractionalWJ);

                    FieldVector<Scalar,dimWorld> vTotal(K);
                    vTotal *= lambda*(pressI - pressJ)/dist;

                    if (hasGravity)
                    {
                        Ki += Kj;
                        Ki *= 0.5;
                        FieldVector<Scalar,dimWorld> gEffect(0);
                        Ki.umv(gravity, gEffect);
                        Scalar factor = fractionalW*(this->diffProblem.wettingPhase.density())
                            + (1 - fractionalW)*(this->diffProblem.nonWettingPhase.density());
                        gEffect *= lambda*factor;
                        vTotal += gEffect;
                    }

                    if (this->diffProblem.capillarity)
                    {
                        // arithmetic average of the permeability
                        K = ((Kni + Knj) *= 0.5);

                        // capillary pressure w.r.t. saturation
                        Scalar pCI = this->diffProblem.materialLaw().pC(satI,globalPos, *eIt, localPos);
                        Scalar pCJ = this->diffProblem.materialLaw().pC(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);

                        // mobility of the nonwetting phase
                        Scalar lambdaN = 0.5*(this->diffProblem.materialLaw().mobN(1 - satI,globalPos, *eIt, localPos)
                                + this->diffProblem.materialLaw().mobN(1 - satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor));

                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= -(pCJ - pCI)/(dist*dist);

                        vTotal -= lambdaN*(K*pCGradient);
                    }

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

                        // compute directed permeability vector Ki.n
                        FieldVector<Scalar,dim> Kni(0);
                        Ki.umv(distVec, Kni);

                        // compute averaged total mobility
                        Scalar lambda = 1.;
                        Scalar fractionalW = 1.;
                        lambda = lambdaI;
                        if (hasGravity) fractionalW = fractionalWI;

                        Scalar g = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);

                        FieldVector<Scalar,dim> vTotal(Kni);
                        vTotal *= lambda*(g-pressI)/dist;
                        if (hasGravity)
                        {
                            FieldVector<Scalar,dimWorld> gEffect(0);
                            Ki.umv(gravity, gEffect);
                            Scalar factor = fractionalW*(this->diffProblem.wettingPhase.density())
                                + (1 - fractionalW)*(this->diffProblem.nonWettingPhase.density());
                            gEffect *= lambda*factor;
                            vTotal += gEffect;
                        }

                        if (this->diffProblem.capillarity)
                        {
                            Scalar satJ = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);

                            // distance vector between barycenters
                            FieldVector<Scalar,dimWorld>
                            distVec = globalPos - globalPosFace;

                            // compute distance between cell centers
                            Scalar dist = distVec.two_norm();

                            // capillary pressure w.r.t. saturation
                            Scalar pCI = this->diffProblem.materialLaw().pC(satI,globalPos, *eIt, localPos);
                            Scalar pCJ = this->diffProblem.materialLaw().pC(satJ,globalPosFace, *eIt, localPosFace);

                            // mobility of the nonwetting phase
                            Scalar lambdaN = 0.5*(this->diffProblem.materialLaw().mobN(1 - satI,globalPos, *eIt, localPos)
                                    + this->diffProblem.materialLaw().mobN(1 - satJ,globalPosFace, *eIt, localPosFace));

                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= -(pCJ - pCI)/(dist*dist);

                            vTotal -= lambdaN*(Kni*pCGradient);
                        }

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
            // end all intersections
            //            std::cout<<"velocity = "<< this->diffProblem.variables.velocity <<std::endl;
            if (dim == 1&& this->diffProblem.capillarity != true)
            {
                Scalar sum = (fabs(this->diffProblem.variables.velocity[globalIdxI][0][0]*faceVol[0])
                              + fabs(this->diffProblem.variables.velocity[globalIdxI][1][0]));
                Scalar diff = fabs(this->diffProblem.variables.velocity[globalIdxI][0][0]*faceVol[0]
                                   - this->diffProblem.variables.velocity[globalIdxI][1][0]*faceVol[1])/sum;
                if (diff> 1e-6&& sum> 1e-9)
                {
                    std::cout << "NOT conservative!!! diff = "<< diff
                              << ", globalIdxI = "<< globalIdxI << std::endl;
                    std::cout << this->diffProblem.variables.velocity[globalIdxI][0][0]*faceVol[0]<< ", "
                              << this->diffProblem.variables.velocity[globalIdxI][1][0]*faceVol[1]<< std::endl;
                }
            }
            if (dim == 2&& this->diffProblem.capillarity != true)
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

        // for test
        //printvector(std::cout, this->diffProblem.variables.velocity, "velocity", "row", 4, 1, 10);

        return;
    }
};
}
#endif
