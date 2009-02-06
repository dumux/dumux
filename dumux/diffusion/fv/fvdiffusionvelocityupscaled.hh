// $Id$
#ifndef DUNE_FVDIFFUSIONVELOCITYUPSCALED_HH
#define DUNE_FVDIFFUSIONVELOCITYUPSCALED_HH

#include "dumux/diffusion/fv/fvdiffusion.hh"
#include <../../../../dune-subgrid/subgrid/subgrid.hh>

namespace Dune
{
/** \todo Please doc me! */

template<class Grid, class Scalar, class VC>
class FVDiffusionVelocityUpscaled: public FVDiffusion<Grid, Scalar, VC>
{

    enum
    {
        dim = Grid::dimension
    };
    enum
    {
        dimWorld = Grid::dimensionworld
    };

typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Grid::template Codim<1>::EntityPointer FacePointer;
    typedef typename GridView::Intersection Intersection;

    typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    FVDiffusionVelocityUpscaled(Grid& grid, FractionalFlowProblem<Grid, Scalar, VC>& problem)
    : FVDiffusion<Grid,Scalar,VC>(grid, problem)
    {}

    void calcTotalVelocity(const Scalar t) const
    {
        int level=this->diffProblem.variables.transLevel;

        if (this->level()<=level)
        {
            DUNE_THROW(NotImplemented,"this->level()=<level");
        }

        const GridView& gridView = this->grid.levelView(level);
        const IndexSet& coarseIndexSet = this->grid.levelIndexSet(level);

        ElementIterator eItCoarseEnd = gridView.template end<0>();
        for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
        {
            int globalIdxICoarse = coarseIndexSet.index(*eItCoarse);

            FieldVector<Scalar, 2*dim> faceAreaCoarse(0);
            BlockVector<FieldVector<Scalar,dim> > fluxCoarse(2*dim);
            fluxCoarse=0;

            HierarchicIterator eItEnd = eItCoarse-> hend(this->level());
            for (HierarchicIterator eIt = eItCoarse->hbegin(this->level()); eIt != eItEnd; ++eIt)
            {
                //only iterat through difflevel!!!
                if (eIt->level() != this->level()) continue;

                //get some cell properties
                GeometryType gt = eIt->geometry().type();
                const LocalPosition&
                localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);
                const GlobalPosition& globalPos = eIt->geometry().global(localPos); //globalPosFace coordinates of cell center
                int globalIdxI = this->elementMapper.map(*eIt); // index of fine-scale cell

                // run through all intersections with neighbors and boundary
                IntersectionIterator isItEnd = gridView.template iend(*eIt);
                for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt!=isItEnd; ++isIt)
                {
                    if (isIt->neighbor())
                    {
                        // neighbor's properties
                        ElementPointer neighborPointer = isIt->outside();
                        ElementPointer neighborPointerCoarse = (*neighborPointer).father();
                        int neighborFatherLevel = neighborPointerCoarse.level();
                        while(level != neighborFatherLevel)
                        {
                            Element& neighborFatherElement = *neighborPointerCoarse;
                            neighborPointerCoarse = neighborFatherElement.father();
                            neighborFatherLevel = neighborPointerCoarse.level();
                        }

                        int globalIdxNeighborCoarse = coarseIndexSet.index(*neighborPointerCoarse);

                        if(globalIdxICoarse != globalIdxNeighborCoarse)
                        {
                            //                            std::cout<<globalIdxICoarse<<globalIdxNeighborCoarse<<std::endl;
                            int faceNumberFine = isIt->numberInSelf();

                            //numbering of the faces is the same for different levels
                            int faceNumberCoarse = faceNumberFine;

                            if (!faceAreaCoarse[faceNumberCoarse])
                            {
                                FacePointer coarseFace = (*eItCoarse).template entity<1>(faceNumberCoarse);
                                faceAreaCoarse[faceNumberCoarse] = (*coarseFace).geometry().volume();
                            }

                            GeometryType faceGT = isIt->intersectionSelfLocal().type();

                            const FieldVector<Scalar,dim-1>& faceLocal
                            = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                            int globalIdxJ = this->elementMapper.map(*neighborPointer);// neigbor's fine-scale cell index

                            GeometryType neighborGT = neighborPointer->geometry().type();

                            const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                            const GlobalPosition&
                            globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor); // globalPosFace coordinate of neighbor's cell center

                            // get pressure and permeability and total mobility in fine-scale element
                            Scalar pressI = this->diffProblem.variables.pressure[globalIdxI];

                            FieldMatrix permeabilityI = this->diffProblem.soil.K(globalPos,*eIt,localPos);

                            Scalar satI = this->diffProblem.variables.saturation[globalIdxICoarse];
                            Scalar lambdaI = this->diffProblem.materialLaw.mobTotal(satI,globalPos,*eIt,localPos);

                            Scalar faceAreaFine;
                            // run through all intersections with neighbors and boundary

                            // get some face properties
                            switch(Grid::dimension)
                            {
                                case 1: faceAreaFine = 1;
                                break;
                                default: faceAreaFine = isIt->intersectionGlobal().volume(); // volume of face
                                break;
                            }

                            FieldVector<Scalar,dim> unitOuterNormal
                            = isIt->unitOuterNormal(faceLocal); // normal vector of unit length

                            // get neighbor pressure and permeability
                            Scalar pressJ = this->diffProblem.variables.pressure[globalIdxJ];

                            FieldMatrix permeabilityJ = this->diffProblem.soil.K(globalPosNeighbor,*neighborPointer,localPosNeighbor);

                            // compute vectorized permeabilities
                            FieldVector<Scalar,dim> normalPermeabilityI(0);
                            FieldVector<Scalar,dim> normalPermeabilityJ(0);
                            permeabilityI.umv(unitOuterNormal, normalPermeabilityI);
                            permeabilityJ.umv(unitOuterNormal, normalPermeabilityJ);
                            // compute permeability normal to intersection and take harmonic mean
                            Scalar normalComponentPermeabilityI = normalPermeabilityI * unitOuterNormal;
                            Scalar normalComponentPermeabilityJ = normalPermeabilityJ * unitOuterNormal;
                            Scalar meanNormalPermeability = 2 * normalComponentPermeabilityI * normalComponentPermeabilityJ / (normalComponentPermeabilityI + normalComponentPermeabilityJ);
                            // compute permeability tangential to intersection and take arithmetic mean
                            FieldVector<Scalar,dim> normalComponentVector = unitOuterNormal;
                            FieldVector<Scalar,dim> tangentialPermeabilityI = normalPermeabilityI - (normalComponentVector *= normalComponentPermeabilityI);
                            normalComponentVector = unitOuterNormal;
                            FieldVector<Scalar,dim> tangentialPermeabilityJ = normalPermeabilityJ - (normalComponentVector *= normalComponentPermeabilityJ);
                            FieldVector<Scalar,dim> meanTangentialPermeability = (tangentialPermeabilityI += tangentialPermeabilityJ);
                            meanTangentialPermeability *= 0.5;
                            FieldVector<Scalar,dim> meanNormalPermeabilityVector = unitOuterNormal;
                            // Build vectorized averaged permeability
                            FieldVector<Scalar,dim> permeability = (meanTangentialPermeability += (meanNormalPermeabilityVector *= meanNormalPermeability));

                            //                            std::cout<<"permeability = "<<permeability<<std::endl;

                            // distance between cell centers
                            FieldVector<Scalar,dimWorld>
                            distVec = globalPos - globalPosNeighbor;
                            Scalar dist = distVec.two_norm();
                            // get averaged total mobility
                            Scalar satJ = this->diffProblem.variables.saturation[globalIdxNeighborCoarse];
                            Scalar lambdaJ = this->diffProblem.materialLaw.mobTotal(satJ,globalPosNeighbor,*neighborPointer,localPosNeighbor);
                            Scalar meanLambda = 0.5 * (lambdaI + lambdaJ);
                            // compute total velocity with Darcy's Law
                            FieldVector<Scalar, dim> velocityFine(permeability);
                            velocityFine *= (meanLambda * (pressI-pressJ) / dist);
                            fluxCoarse[faceNumberCoarse] += (velocityFine*=faceAreaFine);
                        }
                    }

                    // boundary face

                    else
                    {
                        if (isIt->boundary())
                        {
                            GeometryType faceGT = isIt->intersectionSelfLocal().type();

                            const FieldVector<Scalar,dim-1>& faceLocal
                            = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                            const GlobalPosition& globalPosFace= isIt->intersectionGlobal().global(faceLocal); // globalPosFace coordinate of face center

                            int faceNumberFine = isIt->numberInSelf();
                            //numbering of the faces is the same for different levels
                            int faceNumberCoarse = faceNumberFine;

                            if (!faceAreaCoarse[faceNumberCoarse])
                            {
                                FacePointer coarseFace = (*eItCoarse).template entity<1>(faceNumberCoarse);
                                faceAreaCoarse[faceNumberCoarse] = (*coarseFace).geometry().volume();
                            }

                            // get pressure and permeability and total mobility in fine-scale element
                            Scalar pressI = this->diffProblem.variables.pressure[globalIdxI];

                            FieldMatrix permeabilityI = this->diffProblem.soil.K(globalPos,*eIt,localPos);

                            Scalar satI = this->diffProblem.variables.saturation[globalIdxICoarse];
                            Scalar lambdaI = this->diffProblem.materialLaw.mobTotal(satI,globalPos,*eIt,localPos);

                            Scalar faceAreaFine;

                            // get some face properties
                            switch(Grid::dimension)
                            {
                                case 1: faceAreaFine = 1;
                                break;
                                default: faceAreaFine = isIt->intersectionGlobal().volume(); // volume of face
                                break;
                            }

                            FieldVector<Scalar,dim> unitOuterNormal
                            = isIt->unitOuterNormal(faceLocal); // normal vector of unit length

                            const LocalPosition&
                            localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumberFine,1);

                            FieldVector<Scalar,dim> velocityFine(0);

                            //get boundary condition for boundary face center
                            BoundaryConditions::Flags bctype = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);
                            if (bctype == BoundaryConditions::dirichlet)
                            {
                                // distance vector between barycenters
                                FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;
                                Scalar dist = distVec.two_norm();

                                //normalise distVec for multiplication with the permeability
                                distVec /= dist;

                                // compute directed permeability vector permeabilityI.n
                                FieldVector<Scalar,dim> normalPermeabilityI(0);
                                permeabilityI.umv(distVec, normalPermeabilityI);

                                // compute averaged total mobility
                                Scalar meanLambda = lambdaI;

                                Scalar pressBound = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                                //std::cout<<"pressBound = "<<pressBound<<std::endl;
                                //std::cout<<"pressI = "<<pressI<<std::endl;
                                //std::cout<<"normalPermeabilityI = "<<normalPermeabilityI<<std::endl;
                                velocityFine=normalPermeabilityI;
                                velocityFine *= (meanLambda * (pressBound-pressI) / dist);
                                //                                std::cout<<"velocityFine = "<<velocityFine<<std::endl;
                            }
                            else
                            {
                                Scalar neumannFlux = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                                velocityFine=unitOuterNormal;
                                velocityFine *= neumannFlux;
                            }
                            fluxCoarse[faceNumberCoarse] += (velocityFine*=faceAreaFine);
                        }
                    }
                }
                // end intersection traversal
            }// end hierarchic iteration


            // evaluate mean velocities at coarse element edges
            // and write them to the velocity struct.
            for (int i = 0; i<2*dim; i++)
            {
                this->diffProblem.variables.velocity[globalIdxICoarse][i] = fluxCoarse[i]/=faceAreaCoarse[i];
//                                std::cout<<"velocity "<<fluxCoarse[i]<<std::endl;
            }

            // check for conservativity
            if (dim == 2)
            {
                Scalar diff = fabs(this->diffProblem.variables.velocity[globalIdxICoarse][0][0]
                        - this->diffProblem.variables.velocity[globalIdxICoarse][1][0]
                        + this->diffProblem.variables.velocity[globalIdxICoarse][2][1]
                        - this->diffProblem.variables.velocity[globalIdxICoarse][3][1])
                /(fabs(this->diffProblem.variables.velocity[globalIdxICoarse][0][0])
                        + fabs(this->diffProblem.variables.velocity[globalIdxICoarse][1][0])
                        + fabs(this->diffProblem.variables.velocity[globalIdxICoarse][2][1])
                        + fabs(this->diffProblem.variables.velocity[globalIdxICoarse][3][1]));
                if (diff> 1e-6)
                {
                    std::cout << "NOT conservative!!! diff = " << diff << ", globalIdxI = " << globalIdxICoarse << std::endl;
                    std::cout << this->diffProblem.variables.velocity[globalIdxICoarse][0][0] << ", "
                    << this->diffProblem.variables.velocity[globalIdxICoarse][1][0] << ", "
                    << this->diffProblem.variables.velocity[globalIdxICoarse][2][1] << ", "
                    << this->diffProblem.variables.velocity[globalIdxICoarse][3][1] << std::endl;
                }
            }
        } // end grid traversal

        return;
    }// end method totalVelocity

};
}
#endif
