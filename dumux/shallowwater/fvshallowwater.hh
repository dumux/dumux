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

template<class Grid, class Scalar, class VC, class Problem = ShallowProblemBase<Grid, Scalar, VC> > class FVShallowWater :
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
    {   dim = Grid::dimension, oneD = 1, twoD = 2};
    enum
    {   dimWorld = Grid::dimensionworld};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout>
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

    // int reconstructFaceValues()

    void initialize();

    void postProcessUpdate(Scalar t, Scalar dt);

    FVShallowWater(Grid& grid, Problem& problem,
            NumericalFlux<Grid,Scalar>& numFl = *(new HllFlux<Grid,Scalar>)) :
        ShallowWater<Grid, Scalar, VC>(grid, problem),
                elementMapper_(grid.leafView()), numFlux_(numFl),
                gravity_(problem.gravityConstant())
    {
    }

private:
    ElementMapper elementMapper_;
    NumericalFlux<Grid,Scalar>& numFlux_;
    const Scalar& gravity_;
};

//update method computes the update factors for the whole grid in one time step
template<class Grid, class Scalar, class VC, class Problem> int FVShallowWater<Grid, Scalar,
        VC, Problem>::update(const Scalar t, Scalar& dt, SolutionType& updateVec)
{
    // initialize dt very large so model is forced to take the cfl timestep
    dt = 1e100;
    updateVec = 0;

    const GridView& gridView= this->grid.leafView();

    //**************BEGIN GRID TRAVERSAL****************************************************************+

    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        Scalar bottomElevationI=0; //bottom elevation of entitiy
        VelType bottomSlope(0); //slope between two entities (from center to center)
        Scalar waterDepthI=0; //water depth of entity
        Scalar waterLevelI = 0; //water level of entitiy = water depth + bottom elevation
        VelType velocityI(0);
        Scalar maxCflVelocity=0; //decisive velocity component x or y
        Scalar dist=0;
        Scalar froudeNumber = 0;
        VelType divergenceTerm(0); //term computed when applying the divergence form of bed slope term
        VelType summedDivergenceTerm(0); //divergence term summed over faces
        SystemType summedFluxes(0); //summed fluxed over faces
        SystemType sourceTermVector(0);//all sources including rainfall, slope terms and later friction terms
        SystemType flux(0); //resulting flux vector


        //***********GET BASIC VALUES OF ENTITY (geometry, global position, water depth, velociy,...)
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition &localPos =
                Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        GlobalPosition globalPos = eIt->geometry().global(localPos);
        //        std::cout<<"globalPos "<<globalPos<<std::endl;

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)
        *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
        //        std::cout<<volume<<std::endl;

        // cell index
        int globalIdxI = elementMapper_.map(*eIt);
        //        std::cout<<globalIdxI<<std::endl;

        // get bottomElevation
        bottomElevationI = this->problem.variables.bottomElevation[globalIdxI];
        //        std::cout<<"bottomElevationI "<<bottomElevationI<<std::endl;

        // get bottomslope for entity
        // bottomSlope =this->problem.surface.evalBottomSlopes();
        //        std::cout<<"bottom slope"<<bottomSlope<<std::endl;

        // get waterdepth at cell center
        waterDepthI = this->problem.variables.waterDepth[globalIdxI];
        //        std::cout<<"waterDepthI="<<waterDepthI<<std::endl;

        // determine water level of cell
        waterLevelI = waterDepthI + bottomElevationI;

        // get velocity at cell center
        velocityI = this->problem.variables.velocity[globalIdxI];
        //        std::cout<<"velocityI="<<velocityI<<std::endl;

        // Compute Cfl Velocity for Entitiy
        switch (dim)
        {
        case oneD:
            maxCflVelocity = std::max(fabs(velocityI[0])
                    +sqrt(fabs(gravity_*waterDepthI)), maxCflVelocity);
            break;

        case twoD:
            Scalar cflVelocityIX = fabs(velocityI[0])
                    +sqrt(fabs(gravity_*waterDepthI)); //Calculate cfl velocity (convective velocity + wave velocity)
            Scalar cflVelocityIY= fabs(velocityI[1])
                    +sqrt(fabs(gravity_*waterDepthI));

            Scalar cflVelocity = std::max(cflVelocityIX,
                    cflVelocityIY);

            maxCflVelocity = std::max(maxCflVelocity,
                    cflVelocity);
            break;
        }
        // std::cout<<"maxCflVelocity "<<maxCflVelocity<<std::endl;

        // Compute Froude-Number for Entity
        switch (dim)
        {
        case oneD:
            froudeNumber = fabs(velocityI[0])/(sqrt(gravity_*waterDepthI));
            break;
        case twoD:
            Scalar froudeX = fabs(velocityI[0])/(sqrt(gravity_*waterDepthI));
            Scalar froudeY = fabs(velocityI[1])/(sqrt(gravity_*waterDepthI));
            froudeNumber = std::max(froudeX, froudeY);
            break;
        }

        //************BEGIN INTERSECTION TRAVERSAL*************************************************

        Scalar waterDepthJ=0;
        Scalar bottomElevationJ;
        VelType velocityJ(0);
        VelType cflVelocityFace(0);

        IntersectionIterator isItEnd =gridView.template iend(*eIt);
        for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
//            int indexInInside = isIt->indexInInside();//numberInInside

            // get geometry type of face
            Dune::GeometryType gtf = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>& faceLocalPos =
                    Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0, 0);

            // determine volume of face to multiply it with the flux
            Scalar faceVolume = isIt->geometry().volume(); //geometry()
            // std::cout<<faceVolume<<std::endl; //ok

            // center of face inside volume reference element
            const LocalPosition& faceLocalDim =
                    Dune::ReferenceElements<Scalar,dim>::general(gtf).position(isIt->indexInInside(), 1);

            //get normal vector of face
            Dune::FieldVector<Scalar,dimWorld> nVec =
                    isIt->unitOuterNormal(faceLocalPos);
            //            std::cout<<nVec<<std::endl;

            // *********** INTERIOR FACE ********************
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = elementMapper_.map(*neighborPointer);

                Dune::GeometryType nbgt = neighborPointer->geometry().type();
                const LocalPosition& nbLocalPos =
                        Dune::ReferenceElements<Scalar,dim>::general(nbgt).position(0, 0);

                // neighbor cell center in global coordinates
                Dune::FieldVector<Scalar,dimWorld> nbGlobalPos =
                        neighborPointer->geometry().global(nbLocalPos);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos
                        - nbGlobalPos;

                // compute distance between cell centers
                dist = distVec.two_norm();
                //                std::cout<<dist<<std::endl;

                //get bottomElevation
                bottomElevationJ
                        = this->problem.variables.bottomElevation[globalIdxJ];
                // std::cout<<"bottomElevationJ "<<bottomElevationJ<<std::endl;

                // calculate bottom slope
                bottomSlope = (bottomElevationJ - bottomElevationI)/dist;
                //                std::cout<<"bottomslope "<<bottomSlope<<std::endl;

                // get waterdepth at neighbor cell center
                waterDepthJ = this->problem.variables.waterDepth[globalIdxJ];

                // get velocity at neighbor cell center
                velocityJ = this->problem.variables.velocity[globalIdxJ];


                //******INSERT RECONSTRUCTION METHOD reconstructFaceValues() if model
                //has to be of higher order than 1
                //compute variables waterDepthFaceI, waterDepthFaceJ, velocityFaceI, velocityFaceJ
                //and commit it to numerical flux function
                //************************************************************************


                //get flux computed in numerical flux
                flux = (numFlux_(velocityI, velocityJ, waterDepthI,
                        waterDepthJ, nVec));

                flux*=faceVolume;

                // std::cout<<"cell "<<globalIdxI<<"flux_interior of face "<<indexInInside<<" "<<flux<<std::endl;

                // determine bottomelevation at face
                Scalar bottomElevationFace =
                        (bottomElevationI+bottomElevationJ)/2;

                // determine water Depth at face for divergence term purposes
                Scalar waterDepthFaceDivergence = waterLevelI
                        - bottomElevationFace;
                // std::cout<<"waterDepthFaceDivergence "
                //       <<waterDepthFaceDivergence<<std::endl;

                // Compute Divergence Term for Interior Faces
                switch (dim)
                {
                case oneD:
                    divergenceTerm[0] = 0.5*gravity_*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;
                    break;
                case twoD:
                    divergenceTerm[0] = 0.5*gravity_*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;

                    divergenceTerm[1] = 0.5*gravity_ *(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[1]*= nVec[1];
                    divergenceTerm[1]*= faceVolume;
                    break;
                }
                // std::cout<<divergenceTerm<<std::endl;
            }

            //********* BOUNDARY FACE *********************************************

            Scalar waterDepthFace=0;
            VelType velocityFace(0);
            VelType momentumFace(0);
            Scalar contiFlux = 0;
            VelType momentumFlux(0);

            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition faceGlobalPos = isIt->geometry().global(faceLocalPos);
                // std::cout<<"FaceGlobalPos "<<faceGlobalPos<<std::endl;

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos
                        - faceGlobalPos;

                //compute distance between cell centers
                // dist = distVec.two_norm();

                //*************BOUNDARY CONDITIONS *******************

                // get boundary type for conti part
                BoundaryConditions::Flags bctypeConti =
                        this->problem.bctypeConti(faceGlobalPos, *eIt,
                                faceLocalDim, froudeNumber);

                // get boundary type for momentum part
                BoundaryConditions::Flags bctypeMomentum =
                        this->problem.bctypeMomentum(faceGlobalPos, *eIt,
                                faceLocalDim, froudeNumber);

                // possible combinations of boundary conditions for the continuity equation (Conti) and the momentum equations (Momentum)

                if (bctypeConti == BoundaryConditions::dirichlet
                        && bctypeMomentum == BoundaryConditions::dirichlet)
                {
                    waterDepthFace = this->problem.dirichletConti(
                            faceGlobalPos, *eIt, faceLocalDim);

                    momentumFace = this->problem.dirichletMomentum(
                            faceGlobalPos, *eIt, faceLocalDim);

                    velocityFace = momentumFace;
                    velocityFace /= waterDepthFace;

                    switch (dim)
                    {
                    case oneD:

                        contiFlux = velocityFace[0] * waterDepthFace*nVec[0];
                        momentumFlux =(velocityFace[0] * velocityFace[0]
                                * waterDepthFace +0.5 * gravity_
                                * waterDepthFace * waterDepthFace) *nVec[0];
                        break;
                    case twoD:

                        contiFlux = waterDepthFace*velocityFace[0]*nVec[0]
                                +waterDepthFace * velocityFace[1]*nVec[1];

                        momentumFlux[0]= (velocityFace[0] * velocityFace[0]
                                *waterDepthFace +0.5 * gravity_
                                * waterDepthFace *waterDepthFace) * nVec[0]
                                +(waterDepthFace *velocityFace[0]
                                        * velocityFace[1]) * nVec[1];
                        momentumFlux[1]= (waterDepthFace *velocityFace[0]
                                *velocityFace[1]) * nVec[0]+(velocityFace[1]
                                *velocityFace[1] * waterDepthFace + 0.5
                                * gravity_ *waterDepthFace * waterDepthFace)
                                * nVec[1];
                        break;
                    }
                    // std::cout<<"flux dirichlet "<<flux<<std::endl;
                    // assign the conti and momentum fluxes to the overall flux vector
                    flux[0] = contiFlux;

                    for (int i=0; i<dim; i++)
                    {
                        flux[i+1]=momentumFlux[i];
                    }
                    flux*= faceVolume;
                }

                if (bctypeConti == BoundaryConditions::dirichlet
                        && bctypeMomentum == BoundaryConditions::neumann)
                {
                    waterDepthFace = this->problem.dirichletConti(
                            faceGlobalPos, *eIt, faceLocalDim);
                    velocityFace = velocityI;

                    switch (dim)
                    {
                    case oneD:
                        contiFlux = velocityFace[0] * waterDepthFace*nVec[0];
                        momentumFlux =(velocityFace*velocityFace*waterDepthFace
                                +0.5 *gravity_*waterDepthFace*waterDepthFace)
                                *nVec[0];
                        break;
                    case twoD:
                        contiFlux = waterDepthFace*velocityFace[0]*nVec[0]
                                +waterDepthFace*velocityFace[1]*nVec[1];

                        momentumFlux[0]= (velocityFace[0]*velocityFace[0]
                                *waterDepthFace +0.5*gravity_*waterDepthFace
                                *waterDepthFace)*nVec[0] +(waterDepthFace
                                *velocityFace[0]*velocityFace[1]) *nVec[1];
                        momentumFlux[1]= (waterDepthFace *velocityFace[0]
                                *velocityFace[1])*nVec[0]+(velocityFace[1]
                                *velocityFace[1]*waterDepthFace +0.5 * gravity_
                                *waterDepthFace *waterDepthFace)*nVec[1];
                        break;
                    }

                    // Check whether there is a free flow or a no-flow condition imposed
                    momentumFlux = this->problem.neumannMomentum(faceGlobalPos,
                            *eIt, faceLocalDim, waterDepthFace, momentumFlux);

                    // assign the conti and momentum fluxes to the overall flux vector
                    flux[0] = contiFlux;

                    for (int i=0; i<dim; i++)
                    {
                        flux[i+1]=momentumFlux[i];
                    }
                    flux*= faceVolume;
                }

                if (bctypeConti == BoundaryConditions::neumann
                        && bctypeMomentum == BoundaryConditions::dirichlet)
                {
                    momentumFace = this->problem.dirichletMomentum(
                            faceGlobalPos, *eIt, faceLocalDim);

                    waterDepthFace = waterDepthI;

                    velocityFace= momentumFace;
                    velocityFace /= waterDepthFace;

                    switch (dim)
                    {
                    case oneD:
                        contiFlux = waterDepthFace*velocityFace[0]*nVec[0];
                        momentumFlux =(velocityFace[0]*waterDepthFace
                                *velocityFace[0]+0.5 *gravity_ *waterDepthFace
                                *waterDepthFace)*nVec[0];
                        break;
                    case twoD:

                        contiFlux = waterDepthFace*velocityFace[0]*nVec[0]
                                +waterDepthFace*velocityFace[1] *nVec[1];

                        momentumFlux[0]= (velocityFace[0]*waterDepthFace
                                *velocityFace[0] +0.5 *gravity_*waterDepthFace
                                *waterDepthFace)*nVec[0] +(waterDepthFace
                                *velocityFace[0]*velocityFace[1]) *nVec[1];
                        momentumFlux[1]= (waterDepthFace*velocityFace[1]
                                *velocityFace[0] ) *nVec[0]+(velocityFace[1]
                                *waterDepthFace*velocityFace[1] +0.5 * gravity_
                                *waterDepthFace *waterDepthFace)*nVec[1];
                        break;
                    }
                    //assign the conti and momentum fluxes to the overall flux vector
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
                    waterDepthFace = waterDepthI;
                    velocityFace = velocityI;

                    switch (dim)
                    {
                    case oneD:
                        contiFlux = waterDepthFace*velocityFace[0]*nVec[0];
                        momentumFlux =(velocityFace[0]*waterDepthFace
                                *velocityFace[0]+0.5 *gravity_ *waterDepthFace
                                *waterDepthFace)*nVec[0];
                        break;
                    case twoD:
                        contiFlux = waterDepthFace*velocityFace[0]*nVec[0]
                                +waterDepthFace*velocityFace[1] *nVec[1];

                        momentumFlux[0]= (velocityFace[0]*waterDepthFace
                                *velocityFace[0] +0.5 *gravity_*waterDepthFace
                                *waterDepthFace)*nVec[0] +(waterDepthFace
                                *velocityFace[0]*velocityFace[1]) *nVec[1];
                        momentumFlux[1]= (waterDepthFace*velocityFace[1]
                                *velocityFace[0] ) *nVec[0]+(velocityFace[1]
                                *waterDepthFace*velocityFace[1] +0.5 *gravity_
                                *waterDepthFace *waterDepthFace)*nVec[1];
                        break;
                    }
                    // std::cout<<"conti flux "<<contiFlux<<std::endl;

                    // Check whether there is a free flow or a no-flow condition imposed
                    contiFlux = this->problem.neumannConti(faceGlobalPos, *eIt,
                            faceLocalDim, contiFlux);
                    momentumFlux = this->problem.neumannMomentum(faceGlobalPos,
                            *eIt, faceLocalDim, waterDepthFace, momentumFlux);

                    // std::cout<<"momentum flux neumann "<<momentumFlux<<std::endl;

                    // assign the conti and momentum fluxes to the overall flux vector
                    flux[0]= contiFlux;

                    for (int i=0; i<dim; i++)
                    {
                        flux[i+1]=momentumFlux[i];
                    }
                    flux*= faceVolume;
                    // std::cout<<"flux neumann "<<flux<<std::endl;
                }

                // std::cout<<"waterDepthFace " <<waterDepthFace<<std::endl;
                // std::cout<<"velocityFace " <<velocityFace<<std::endl;
                // std::cout<<flux<<std::endl;


                // Compute Divergence Term for Boundary Faces
                // assumption: bottomelevation at face is same as of cell
                Scalar bottomElevationFace = bottomElevationI;

                // determine water Depth at face for divergence term purposes
                Scalar waterDepthFaceDivergence = waterLevelI
                        - bottomElevationFace;
                // std::cout<<"waterDepthFaceDivergence "
                //       <<waterDepthFaceDivergence<<std::endl;

                switch (dim)
                {
                case oneD:
                    divergenceTerm[0] = 0.5*gravity_*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;
                    break;
                case twoD:
                    divergenceTerm[0] = 0.5*gravity_*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[0]*= nVec[0];
                    divergenceTerm[0]*= faceVolume;

                    divergenceTerm[1] = 0.5*gravity_*(waterDepthFaceDivergence
                            *waterDepthFaceDivergence);
                    divergenceTerm[1]*= nVec[1];
                    divergenceTerm[1]*= faceVolume;
                    break;
                }

                // Determine cflVelocity at faces and store the maximum
                switch (dim)
                {
                case oneD:
                    maxCflVelocity = std::max(maxCflVelocity, fabs(velocityFace[0])
                            +sqrt(fabs(gravity_ *waterDepthFace)));
                    break;

                case twoD:
                    Scalar cflVelocityFaceX= fabs(velocityFace[0])
                            +sqrt(fabs(gravity_ *waterDepthFace)); //Calculate cfl velocity (convective velocity + wave velocity)
                    Scalar cflVelocityFaceY= fabs(velocityFace[1])
                            +sqrt(fabs(gravity_ *waterDepthFace));
                    Scalar cflVelocity = std::max(cflVelocityFaceX, cflVelocityFaceY);

                    maxCflVelocity = std::max(maxCflVelocity,
                            cflVelocity);
                    break;
                }
            }

            // sum the fluxes over the faces
            summedFluxes += flux;
            // sum divergence term over faces
            summedDivergenceTerm+=divergenceTerm;
        }
        //*************END INTERSECTION TRAVERSAL******************************************l

        // std::cout<<"maxCflVelocityFace "<< maxCflVelocityFace
        //                      <<std::endl;
        //   std::cout<<"summedfluxes of cell_"<<globalIdxI<<" = "<<summedFluxes
        //     <<std::endl;
        // std::cout<<"summedDivergenceTerm "<<globalIdxI<<" "
        //    <<summedDivergenceTerm <<std::endl;


        //Get Source Term (rainfall) from varibale class and combine it with divergence term (bed slope)
        //rainfall and bedslope are the two components of teh source term vector
        Scalar sourceTerm = this->problem.setSource(globalPos, *eIt, localPos);

        switch (dim)
        {
        case oneD:
            sourceTermVector[0]=sourceTerm;
            sourceTermVector[1]=summedDivergenceTerm[0];
            break;
        case twoD:
            sourceTermVector[0]=sourceTerm;
            sourceTermVector[1]=summedDivergenceTerm[0];
            sourceTermVector[2]=summedDivergenceTerm[1];
            break;
        }

        // std::cout<<"Source Term Vector of cell"<<globalIdxI<<"= "
        //     <<sourceTermVector<<std::endl;


        //*********** UPDATE CELL VALUES *******************************************************************

        updateVec[globalIdxI] -= summedFluxes;
        updateVec[globalIdxI] += sourceTermVector;
        updateVec[globalIdxI] /= volume;

        // std::cout<<"update for cell"<<globalIdxI<<" = "<<updateVec[globalIdxI]
        //   <<std::endl;

        //*****determine dt with cfl criterium considering the velocities of entities and faces ***********
        Scalar cflTimeStep = 0;

        cflTimeStep = dist/maxCflVelocity;
        dt = std::min(dt, cflTimeStep);

        /* std::cout<<"finalCflVelocity "<< finalCflVelocity<<std::endl;
         std::cout<<"cflTimeStep "<<cflTimeStep<<std::endl;
         std::cout<<"dist "<<dist<<std::endl;
         std::cout<<"dt "<<dt<<std::endl;
         */
    }

    // std::cout<<"dt "<<dt<<std::endl;
    //**************END GRID TRAVERSAL ****************************************************************************************************+

    return 0;
}

//initialize method gets initial values from variableclass
template<class Grid, class Scalar, class VC, class Problem> void FVShallowWater<Grid, Scalar,
        VC, Problem>::initialize()
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

        int globalIdx = elementMapper_.map(*eIt);

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

// postProcessUpdate changes the results of the update method which are
// conserved variables (h, hu, hv) in primary variables (h, u ,v)

template<class Grid, class Scalar, class VC, class Problem> void FVShallowWater<Grid, Scalar,
        VC, Problem>::postProcessUpdate(Scalar t, Scalar dt)
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

        // std::cout<<"global Solution of cell"<<i<<"= "
        //           <<this->problem.variables.globalSolution[i]<<std::endl;

        // std::cout<<"waterDepth = "<<this->problem.variables.waterDepth[i]
        //             <<std::endl;
        // std::cout<<"velX = "<<this->problem.variables.velocity[i][0]<<std::endl;
        // std::cout<<"velY = "<<this->problem.variables.velocity[i][1]<<std::endl;

    }
    return;
}

}
#endif

// end namespace


