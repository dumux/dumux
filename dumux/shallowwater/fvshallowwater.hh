// $Id$

#ifndef DUNE_FVSHALLOWWATER_HH
#define DUNE_FVSHALLOWWATER_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/shallowwater/shallowwater.hh"
#include "dumux/shallowwater/shallownumericalflux.hh"

namespace Dune
{
//! \ingroup transport

template<class Grid, class Scalar, class VC> class FVShallowWater :
    public ShallowWater< Grid, Scalar, VC>
{
    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
    {   dim = Grid::dimension};
    enum
    {   dimWorld = Grid::dimensionworld};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout>
            ElementMapper;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar,dim> VelType;
    typedef Dune::FieldVector<Scalar,dim+1> SystemType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> > SolutionType;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> >
            RepresentationType;

    int update(const Scalar t, Scalar& dt, SolutionType& updateVec);

    void initialize();

    void postProcessUpdate(Scalar t, Scalar dt);

    FVShallowWater(Grid& grid, ShallowProblemBase<Grid, Scalar, VC>& problem,
            NumericalFlux<Grid,Scalar>& numFl = *(new HllFlux2d<Grid,Scalar>)) :
        ShallowWater<Grid, Scalar, VC>(grid, problem), elementMapper(grid,
                grid.leafIndexSet()), numFlux(numFl)
    {
    }

private:
    ElementMapper elementMapper;
    //const IndexSet& indexSet;
    NumericalFlux<Grid,Scalar>& numFlux;

};

template<class Grid, class Scalar, class VC> int FVShallowWater<Grid, Scalar,
        VC>::update(const Scalar t, Scalar& dt, SolutionType& updateVec)
{
    // initialize dt very large, why?
    dt = 1e100;

    updateVec = 0;
    // initialize and declare variables
    const GridView& gridView= this->grid.leafView();

    // compute update vector
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {

        Scalar waterDepthI=0;
        Scalar waterLevelI = 0;
        Scalar bottomElevationI=0;
        Scalar gravity=9.81;
        Scalar dist=0;
        Scalar froudeNumber = 0;

        VelType velocityI(0);
        VelType bottomSlope(0);
        VelType bottomSlopeVector(0);
        VelType divergenceTerm(0);
        VelType summedDivergenceTerm(0);

        SystemType summedFluxes(0);
        SystemType sourceTermVector(0);//all sources including rainfall, slope terms and later on friction terms
        SystemType flux(0);

        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition &localPos =
                Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        GlobalPosition globalPos = eIt->geometry().global(localPos);
        //std::cout<<"globalPos "<<globalPos<<std::endl;

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)
                *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
        //  std::cout<<volume<<std::endl;

        // cell index
        int globalIdxI = elementMapper.map(*eIt);
        //std::cout<<globalIdxI<<std::endl;

        //get bottomElevation
        bottomElevationI = this->problem.variables.bottomElevation[globalIdxI];
        //std::cout<<"bottomElevationI "<<bottomElevationI<<std::endl;

        // get waterdepth at cell center
        waterDepthI = this->problem.variables.waterDepth[globalIdxI];
        //std::cout<<"waterDepthI="<<waterDepthI<<std::endl;

        // determine water level of cell
        waterLevelI = waterDepthI + bottomElevationI;

        // get velocity at cell center
        velocityI = this->problem.variables.velocity[globalIdxI];
        std::cout<<"velocityI="<<velocityI<<std::endl;

        //determine Froude Number of cell
        //  froudeNumber = fabs(velocityI);
        //froudeNumber /= sqrt(gravity*waterDepthI);
        // std::cout<<"froudeNumber of Cell"<<globalIdxI<<" "<<froudeNumber<<std::endl;

        //get bottomslopevector for entity
        // bottomSlope =this->problem.surface.evalBottomSlopes();
        //std::cout<<"bottom slope"<<bottomSlope<<std::endl;


        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd =gridView.template iend(*eIt);
        for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {

            // local number of facet
            //    int numberInSelf = isIt->numberInSelf();//numberInInside

            // get geometry type of face
            Dune::GeometryType gtf = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>& faceLocalPos =
                    Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0, 0);

            //determine volume of face to multiply it with the flux
            Scalar faceVolume = isIt->intersectionGlobal().volume(); //geometry() 
            //  std::cout<<faceVolume<<std::endl; //ok

            // center of face inside volume reference element
            const LocalPosition& faceLocalDim =
                    Dune::ReferenceElements<Scalar,dim>::general(gtf).position(isIt->numberInSelf(), 1);

            //get normal vector of face
            Dune::FieldVector<Scalar,dimWorld> nVec =
                    isIt->unitOuterNormal(faceLocalPos);
            // std::cout<<nVec<<std::endl;

            Scalar waterDepthJ=0;
            Scalar waterDepthFace=0;
            Scalar waterDepthFaceDivergence=0;
            Scalar bottomSlope(0);
            Scalar bottomElevationJ;
            Scalar bottomElevationFace=0;

            VelType velocityJ(0);
            VelType velocityFace(0);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = elementMapper.map(*neighborPointer);

                Dune::GeometryType nbgt = neighborPointer->geometry().type();
                const LocalPosition& nbLocalPos =
                        Dune::ReferenceElements<Scalar,dim>::general(nbgt).position(0, 0);

                // neighbor cell center in global coordinates
                Dune::FieldVector<Scalar,dimWorld> nbGlobalPos =
                        neighborPointer->geometry().global(nbLocalPos);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos
                        - nbGlobalPos;

                //compute distance between cell centers
                dist = distVec.two_norm();
                //std::cout<<dist<<std::endl;    

                //get bottomElevation
                bottomElevationJ
                        = this->problem.variables.bottomElevation[globalIdxJ];
                // std::cout<<"bottomElevationJ "<<bottomElevationJ<<std::endl;

                //calculate bottom slope 
                bottomSlope = (bottomElevationJ - bottomElevationI)/dist;
                //std::cout<<"bottomslope "<<bottomSlope<<std::endl;

                // get waterdepth at neighbor cell center
                waterDepthJ = this->problem.variables.waterDepth[globalIdxJ];

                // get velocity at neighbor cell center
                velocityJ = this->problem.variables.velocity[globalIdxJ];

                //fluxVector for a given scheme (Upwind or different will be catched from numericalflux)
                flux = (numFlux(velocityI, velocityJ, waterDepthI, waterDepthJ,
                        nVec));

                flux*=faceVolume;

                //    std::cout<<"cell "<<globalIdxI<<"flux_interior of face "
                //   <<numberInSelf<<"=" <<flux<<std::endl;

                //******************************************************************************************

                //these steps are only relevant if divergence term for bedslope is considered 

                //determine bottomelevation at face
                Scalar bottomElevationFace =
                        (bottomElevationI+bottomElevationJ)/2;

                //determine water Depth at face for divergence term purposes
                Scalar waterDepthFaceDivergence = waterLevelI
                        - bottomElevationFace;
                // std::cout<<"waterDepthFaceDivergence "
                //       <<waterDepthFaceDivergence<<std::endl;

                //determine divergence Term
                switch (dim)
                {
                case 1:
                    divergenceTerm[0] = 0.5*gravity*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;
                    break;
                case 2:
                    divergenceTerm[0] = 0.5*gravity*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;

                    divergenceTerm[1] = 0.5*gravity *(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[1]*= nVec[1];
                    divergenceTerm[1]*= faceVolume;
                    break;
                }

                //std::cout<<divergenceTerm<<std::endl;

                //******************************************************************************************

            }

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition faceGlobalPos = isIt->intersectionGlobal().global(faceLocalPos);
                // std::cout<<"FaceGlobalPos "<<faceGlobalPos<<std::endl;

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos
                        - faceGlobalPos;

                //compute distance between cell centers
                dist = distVec.two_norm();

                //get boundary type for conti part*************************************************+++
                BoundaryConditions::Flags bctypeConti =
                        this->problem.bctypeConti(faceGlobalPos, *eIt,
                                faceLocalDim);

                //get boundary type for momentum part*************************************************+++
                BoundaryConditions::Flags bctypeMomentum =
                        this->problem.bctypeMomentum(faceGlobalPos, *eIt,
                                faceLocalDim);

                // std::cout<<"bctypeConti"<<bctypeConti <<std::endl;
                //  std::cout<<"bctypeMomentum "<<bctypeMomentum<<std::endl;

                // if waterdepth and momentum are given, velocity can be determined

                if (bctypeConti == BoundaryConditions::dirichlet
                        && bctypeMomentum == BoundaryConditions::dirichlet)
                {
                    waterDepthJ = this->problem.dirichletConti(faceGlobalPos,
                            *eIt, faceLocalDim);

                    VelType momentumJ(0);

                    momentumJ = this->problem.dirichletMomentum(faceGlobalPos,
                            *eIt, faceLocalDim);

                    velocityJ = momentumJ;
                    velocityJ /= waterDepthJ;

                    switch (dim)
                    {
                    case 1:

                        flux[0] = momentumJ[0]*nVec[0];

                        flux[1] = velocityJ[0]*momentumJ[0]+0.5 *9.81
                                *waterDepthJ *waterDepthJ;
                        break;
                    case 2:

                        flux[0] = momentumJ[0]*nVec[0]+momentumJ[1]*nVec[1];

                        flux[1]= (velocityJ[0]*momentumJ[0] +0.5*9.81
                                *waterDepthJ*waterDepthJ)*nVec[0]
                                +(momentumJ[0]*velocityJ[1]) *nVec[1];
                        flux[2]= (momentumJ[1] *velocityJ[0] )*nVec[0]
                                +(velocityJ[1] *momentumJ[0] +0.5*9.8
                                        *waterDepthJ *waterDepthJ)*nVec[1];
                        break;
                    }

                    flux*= faceVolume;
                }

                if (bctypeConti == BoundaryConditions::dirichlet
                        && bctypeMomentum == BoundaryConditions::neumann)
                {
                    waterDepthJ = this->problem.dirichletConti(faceGlobalPos,
                            *eIt, faceLocalDim);
                    velocityJ = velocityI;

                    Scalar contiFlux(0);
                    VelType momentumFlux(0);

                    switch (dim)
                    {
                    case 1:

                        contiFlux = velocityJ[0] * waterDepthJ*nVec[0];
                        momentumFlux =(velocityJ*velocityJ*waterDepthJ+0.5
                                *9.81*waterDepthJ*waterDepthJ)*nVec[0];
                        break;
                    case 2:

                        contiFlux = waterDepthJ*velocityJ[0]*nVec[0]
                                +waterDepthJ*velocityJ[1]*nVec[1];

                        momentumFlux[0]= (velocityJ[0]*velocityJ[0]*waterDepthJ
                                +0.5*9.81*waterDepthJ*waterDepthJ)*nVec[0]
                                +(waterDepthJ *velocityJ[0]*velocityJ[1])
                                        *nVec[1];
                        momentumFlux[1]= (waterDepthJ *velocityJ[0]
                                *velocityJ[1])*nVec[0]+(velocityJ[1]
                                *velocityJ[1]*waterDepthJ +0.5*9.8*waterDepthJ
                                *waterDepthJ)*nVec[1];
                        break;
                    }

                    /*   contiFlux = this->problem.neumannConti(faceGlobalPos, *eIt,
                     faceLocalDim);
                     */
                    momentumFlux = this->problem.neumannMomentum(globalPos,
                            *eIt, faceLocalDim, momentumFlux);

                    flux[0]=contiFlux;

                    for (int i=0; i<dim; i++)
                    {
                        flux[i+1]=momentumFlux[i];
                    }

                    flux*= faceVolume;
                }

                if (bctypeConti == BoundaryConditions::neumann
                        && bctypeMomentum == BoundaryConditions::dirichlet)
                {
                    Scalar contiFlux(0);
                    VelType momentumFlux(0);

                    VelType contiFaceFlux(0);

                    contiFaceFlux = this->problem.neumannConti(faceGlobalPos,
                            *eIt, faceLocalDim);

                    waterDepthJ = waterDepthI;

                    velocityJ= contiFaceFlux;
                    velocityJ /= waterDepthJ;

                    switch (dim)
                    {
                    case 1:

                        contiFlux = contiFaceFlux[0]*nVec[0];
                        momentumFlux =(velocityJ*contiFaceFlux+0.5 *9.81
                                *waterDepthJ *waterDepthJ)*nVec[0];
                        break;
                    case 2:

                        contiFlux = contiFaceFlux[0]*nVec[0]+contiFaceFlux[1]
                                *nVec[1];

                        momentumFlux[0]= (velocityJ[0]*contiFaceFlux[0] +0.5
                                *9.81*waterDepthJ*waterDepthJ)*nVec[0]
                                +(contiFaceFlux[0]*velocityJ[1]) *nVec[1];
                        momentumFlux[1]= (contiFaceFlux[1] *velocityJ[0] )
                                *nVec[0]+(velocityJ[1] *contiFaceFlux[1] +0.5
                                *9.8*waterDepthJ *waterDepthJ)*nVec[1];
                        break;
                    }
                    flux[0]= contiFlux;

                    for (int i=0; i<dim; i++)
                    {
                        flux[i+1]=momentumFlux[i];
                    }

                    flux*= faceVolume;

                }

                if (bctypeConti == BoundaryConditions::neumann
                        && bctypeMomentum == BoundaryConditions::neumann)
                {
                    waterDepthJ = waterDepthI;
                    velocityJ = velocityI;

                    VelType unityVector(0);
                    unityVector[0]= nVec[0]*nVec[0];
                    unityVector[1]= nVec[1]*nVec[1];

                    flux = numFlux(velocityI, velocityJ, waterDepthI,
                            waterDepthJ, nVec);

                    VelType contiFlux(flux[0]);
                    contiFlux = this->problem.neumannConti(globalPos, *eIt,
                            faceLocalDim, contiFlux);

                    flux[0] = contiFlux*unityVector;

                    VelType momentumFlux(0);
                    for (int i=0; i<dim; i++)
                    {
                        momentumFlux[i]=flux[i+1];
                    }

                    momentumFlux = this->problem.neumannMomentum(globalPos,
                            *eIt, faceLocalDim, momentumFlux);

                    for (int i=0; i<dim; i++)
                    {
                        flux[i+1]=momentumFlux[i];
                    }

                    flux*= faceVolume;

                }
                // std::cout<<"waterDepthFace " <<waterDepthFace<<std::endl;
                // std::cout<<"velocityFace " <<velocityFace<<std::endl;


                //std::cout<<flux<<std::endl;

                //determine divergence Term
                //assumption: bottomelevation at face is same as of cell
                Scalar bottomElevationFace = bottomElevationI;

                //determine water Depth at face for divergence term purposes
                Scalar waterDepthFaceDivergence = waterLevelI
                        - bottomElevationFace;
                // std::cout<<"waterDepthFaceDivergence "
                //       <<waterDepthFaceDivergence<<std::endl;

                //determine divergence Term
                switch (dim)
                {
                case 1:
                    divergenceTerm[0] = 0.5*gravity*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;
                    break;
                case 2:
                    divergenceTerm[0] = 0.5*gravity*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;

                    divergenceTerm[1] = 0.5*gravity*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[1]*= nVec[1];
                    divergenceTerm[1]*= faceVolume;
                    break;
                }
            }

            //  std::cout<<flux<<std::endl;

            summedFluxes += flux;
            summedDivergenceTerm+=divergenceTerm;

        }

        // std::cout<<"summedDivergenceTerm "<<globalIdxI<<" "
        //     <<summedDivergenceTerm <<std::endl;


        //     std::cout<<"summedfluxes of cell_"<<globalIdxI<<" = "<<summedFluxes
        //    <<std::endl;

        //Quelle einbinden und in Systemvektor konvertieren
        Scalar sourceTerm = this->problem.setSource(globalPos, *eIt, localPos);

        switch (dim)
        {
        case 1:
            sourceTermVector[0]=sourceTerm;
            sourceTermVector[1]=summedDivergenceTerm[0];
            break;
        case 2:
            sourceTermVector[0]=sourceTerm;
            sourceTermVector[1]=summedDivergenceTerm[0];
            sourceTermVector[2]=summedDivergenceTerm[1];
            break;
        }

        //   std::cout<<"Source Term Vector of cell"<<globalIdxI<<"= "
        //     <<sourceTermVector<<std::endl;

        //Update of cell values************************************************************** 

        updateVec[globalIdxI] -= summedFluxes;
        updateVec[globalIdxI] += sourceTermVector;
        updateVec[globalIdxI] /= volume;

        //       std::cout<<"update for cell"<<globalIdxI<<" = "<<updateVec[globalIdxI]
        //            <<std::endl;

        //cfl-criterium**********************************************************************

        VelType cflVelocity(0);
        Scalar finalCflVelocity = 0;
        Scalar cflTimeStep = 0;

        switch (dim)
        {
        case 1:
            cflVelocity[0]= fabs(velocityI[0])+sqrt(fabs(gravity*waterDepthI));
            finalCflVelocity = cflVelocity[0];
            break;

        case 2:
            cflVelocity[0]= fabs(velocityI[0])+sqrt(fabs(gravity*waterDepthI)); //Calculate cfl velocity (convective velocity + wave velocity)
            cflVelocity[1]= fabs(velocityI[1])+sqrt(fabs(gravity*waterDepthI));

            finalCflVelocity = std::max(cflVelocity[0], cflVelocity[1]);
            break;
        }

        cflTimeStep = dist/finalCflVelocity;
        dt = std::min(dt, cflTimeStep); //calculate time step with cfl criterium

        // std::cout<<"cflvelocity "<< cflVelocity<<std::endl;
        // std::cout<<"cflTimeStep "<< cflTimeStep<<std::endl;
        // std::cout<<"dt "<<dt<<std::endl;
        // std::cout<<"dist "<< dist<<std::endl;

    }

    // end grid traversal

    return 0;
}

template<class Grid, class Scalar, class VC> void FVShallowWater<Grid, Scalar,
        VC>::initialize()
{

    const GridView& gridView= this->grid.leafView();

    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {

        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition &localPos =
                Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        int globalIdx = elementMapper.map(*eIt);

        Scalar initialWaterDepth=this->problem.setInitialWaterDepth(globalPos,
                *eIt, localPos);
        VelType initialVelocity=this->problem.setInitialVelocity(globalPos,
                *eIt, localPos);
        Scalar initialBottomElevation =
                this->problem.surface.evalBottomElevation(globalPos);

        Scalar initialWaterLevel = initialWaterDepth+initialBottomElevation;

        // initialize cell values
        this->problem.variables.waterDepth[globalIdx] = initialWaterDepth;
        this->problem.variables.velocity[globalIdx] = initialVelocity;
        this->problem.variables.globalSolution[globalIdx][0]
                = initialWaterDepth;
        this->problem.variables.globalSolution[globalIdx][1]
                = initialWaterDepth*initialVelocity[0];
        this->problem.variables.globalSolution[globalIdx][dim]
                = initialWaterDepth*initialVelocity[dim-1];
        this->problem.variables.bottomElevation[globalIdx]
                = initialBottomElevation;
        this->problem.variables.waterLevel[globalIdx] = initialWaterLevel;

    }
    return;
}

template<class Grid, class Scalar, class VC> void FVShallowWater<Grid, Scalar,
        VC>::postProcessUpdate(Scalar t, Scalar dt)
{
    for (int i=0; i<this->problem.variables.size; i++)
    {

        if (this->problem.variables.globalSolution[i][0]> 0)
        {

            this->problem.variables.waterDepth[i]
                    =this->problem.variables.globalSolution[i][0];
            this->problem.variables.velocity[i][0]
                    =this->problem.variables.globalSolution[i][1]
                            /this->problem.variables.waterDepth[i];
            this->problem.variables.velocity[i][dim-1]
                    =this->problem.variables.globalSolution[i][dim]
                            /this->problem.variables.waterDepth[i];
            this->problem.variables.waterLevel[i]
                    =this->problem.variables.waterDepth[i]
                            +this->problem.variables.bottomElevation[i];
        }
        else
        {

            this->problem.variables.waterDepth[i]=0;
            this->problem.variables.velocity[i][0]=0;
            this->problem.variables.velocity[i][dim-1]=0;

        }

        //      std::cout<<"global Solution of cell"<<i<<"= "
        //           <<this->problem.variables.globalSolution[i]<<std::endl;

        //      std::cout<<"waterDepth = "<<this->problem.variables.waterDepth[i]
        //             <<std::endl;
        //     std::cout<<"velX = "<<this->problem.variables.velocity[i][0]<<std::endl;
        // std::cout<<"velY = "<<this->problem.variables.velocity[i][1]<<std::endl;

    }
    return;
}

}
#endif

// end namespace


