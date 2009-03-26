// $Id$

#ifndef DUNE_FVSHALLOWWATER_HH
#define DUNE_FVSHALLOWWATER_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/shallowwater/shallowwater.hh"
#include "dumux/shallowwater/shallownumericalflux.hh"
#include "dumux/shallowwater/shallowproblemplain.hh"
#include "dumux/shallowwater/shallowvariableclass.hh"

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

    typedef Dune::FieldVector<Scalar,dim> SlopeType;
    typedef Dune::FieldVector<Scalar,dim+1> SystemType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> > SolutionType;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> >
            RepresentationType;

    int update(const Scalar t, Scalar& dt, SolutionType& updateVec);

    void initialize();

    void postProcessUpdate(Scalar t, Scalar dt);

    FVShallowWater(Grid& grid, ShallowProblemBase<Grid, Scalar, VC>& problem,
            NumericalFlux<Grid,Scalar>& numFl = *(new HllFlux<Grid,Scalar>)) :
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
        Scalar dist=0;
        Scalar waterDepthI=0;
        Scalar waterDepthJ=0;
        Scalar waterDepthFace=0;
        //Scalar waterDepthFaceI;
        //Scalar waterDepthFaceJ;
        SlopeType velocityI(0);
        SlopeType velocityJ(0);
        SlopeType velocityFace(0);
        SlopeType bottomSlope(0); //real bottomslope
        SystemType sourceTermVector(0);//all sources including rainfall, slope terms and later on friction terms
        SystemType flux(0);
        SlopeType divergenceTerm(0);
        SystemType summedFluxes(0);
        Scalar gravity=9.81;
        Scalar cflTimeStep=0;

        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition &localPos =
                Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)
                *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        //  std::cout<<volume<<std::endl;

        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        //get bottomslopevector for entity
        bottomSlope =this->problem.surface.evalBottomSlopes();
        //std::cout<<"bottom slope"<<bottomSlope<<std::endl;

        // get waterdepth at cell center
        waterDepthI = this->problem.variables.waterDepth[globalIdxI];
        //std::cout<<"waterDepthI="<<waterDepthI<<std::endl;

        // get velocity at cell center
        velocityI = this->problem.variables.velocity[globalIdxI];
        //std::cout<<"velocityI="<<velocityI<<std::endl;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd =gridView.template iend(*eIt);
        for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {

            // local number of facet
            int numberInSelf = isIt->numberInSelf();

            // get geometry type of face
            Dune::GeometryType gtf = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>& faceLocalPos =
                    Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0, 0);

            //determine volume of face to multiply it with the flux
      //      Scalar faceVolume = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).volume();
           
                    Scalar faceVolume = isIt->intersectionGlobal().volume();   
                   //  std::cout<<faceVolume<<std::endl; //ok

            // center of face inside volume reference element
            const LocalPosition& faceLocalDim =
                    Dune::ReferenceElements<Scalar,dim>::general(gtf).position(isIt->numberInSelf(), 1);

            //get normal vector of face
            Dune::FieldVector<Scalar,dimWorld> nVec =
                    isIt->unitOuterNormal(faceLocalPos);

            // std::cout<<nVec<<std::endl; // ok


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

                // std::cout<<dist<<std::endl;

                // get waterdepth at neighbor cell center
                waterDepthJ = this->problem.variables.waterDepth[globalIdxJ];

                // get velocity at neighbor cell center
                velocityJ = this->problem.variables.velocity[globalIdxJ];

                //fluxVector for a given scheme (Upwind or different will be catched from numericalflux)

                flux = (numFlux(velocityI, velocityJ, waterDepthI, waterDepthJ,
                        nVec));
                flux*=faceVolume;

                //   std::cout<<"cell "<<globalIdxI<<"flux_interior of face "
                //    <<numberInSelf<<"=" <<flux<<std::endl;

            }
            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition faceGlobalPos = isIt->intersectionGlobal().global(faceLocalPos);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos
                        - faceGlobalPos;

                //compute distance between cell centers
                dist = distVec.two_norm();

                //get boundary type
                BoundaryConditions::Flags bctype = this->problem.bctype(
                        faceGlobalPos, *eIt, faceLocalDim);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    // get waterdepth at boundary
                    //  waterDepthFace
                    //      = this->problem.dirichletWaterDepth(faceGlobalPos);

                    velocityFace
                            = this->problem.dirichletVelocity(faceGlobalPos);

                    waterDepthFace=waterDepthI;

                    flux = numFlux(velocityI, velocityFace, waterDepthI,
                            waterDepthFace, nVec);
                    flux*=faceVolume;

                    //    std::cout<<"flux_dirichlet of face "<<numberInSelf<<"="
                    //  <<flux<<std::endl;

                }
                if (bctype == BoundaryConditions::neumann) //free flow through boundary

                {
                    waterDepthFace = waterDepthI;
                    velocityFace = velocityI;

                    SystemType helpFlux = numFlux(velocityI, velocityFace,
                            waterDepthI, waterDepthFace, nVec);

                    flux = this->problem.neumannFlux(faceGlobalPos, helpFlux);
                    flux*=faceVolume;

                    //     std::cout<<"flux_neumann of face "<<numberInSelf<<"="<<flux
                    //       <<std::endl;
                }

            }
            //   std::cout<<flux<<std::endl;
            summedFluxes += flux;
        }

        //  std::cout<<"summedfluxes of cell_"<<globalIdxI<<" = "<<summedFluxes
        //        <<std::endl;

        //Quelle einbinden und in Systemvektor konvertieren
        Scalar sourceTerm = 0;//this->problem.setSource(globalPos, *eIt, localPos);

        //  switch (dim)
        //  {
        //  case 1:
        sourceTermVector[0]=sourceTerm;
        sourceTermVector[1]=9.81;
        sourceTermVector[1]*=bottomSlope[0]*(-1);
        sourceTermVector[1]*=(waterDepthI);

        //  break;

        /*  case 2:
         sourceTermVector[0]=sourceTerm;
         sourceTermVector[1]=9.81;
         sourceTermVector[1]*=bottomSlope[0]*(-1);
         sourceTermVector[1]*=(waterDepthI);

         sourceTermVector[2]=9.81;
         sourceTermVector[2]*=bottomSlope[1]*(-1);
         sourceTermVector[2]*=(waterDepthI);

         break;
         }*/

        //      std::cout<<"Source Term Vector of cell"<<globalIdxI<<"= "
        //          <<sourceTermVector<<std::endl;

        updateVec[globalIdxI]+=summedFluxes;
        updateVec[globalIdxI]+= sourceTermVector;
        updateVec[globalIdxI]/=(volume);

        //   std::cout<<"update for cell"<<globalIdxI<<" = "<<updateVec[globalIdxI]
        //          <<std::endl;

        //cfl-criterium

        Scalar cflVelocity;

        cflVelocity = fabs(velocityI[0])+sqrt(fabs(gravity*waterDepthI)); //Calculate cfl velocity (convective velocity + wave velocity)
        //   cflVelocity +=sqrt(gravity*waterDepthI);

 //       if (cflVelocity > 0)
   //     {
            cflTimeStep = dist/cflVelocity;            
            dt = std::min(dt, cflTimeStep); //calculate time step with cfl criterium

           
     //   }
 //       else
  //      {
    //        dt = 0.001;
      //  }
         std::cout<<"cflvelocity "<< cflVelocity<<std::endl;
   //      std::cout<<"cflTimeStep "<< cflTimeStep<<std::endl;
     //    std::cout<<"dt "<<dt<<std::endl;
        //  std::cout<<"dist "<< dist<<std::endl;
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
        SlopeType initialVelocity=this->problem.setInitialVelocity(globalPos,
                *eIt, localPos);

        // initialize cell values
        this->problem.variables.waterDepth[globalIdx] = initialWaterDepth;
        this->problem.variables.velocity[globalIdx] = initialVelocity;
        this->problem.variables.globalSolution[globalIdx][0]
                = initialWaterDepth;
        this->problem.variables.globalSolution[globalIdx][1]
                = initialWaterDepth*initialVelocity[0];
        // this->problem.variables.globalSolution[globalIdx][2]
        //     = initialWaterDepth*initialVelocity[1];
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
            // this->problem.variables.velocity[i][1]
            //       =this->problem.variables.globalSolution[i][2]
            //           /this->problem.variables.waterDepth[i];
        }
        else
        {

            this->problem.variables.waterDepth[i]=0;
            this->problem.variables.velocity[i][0]=0;
            // this->problem.variables.velocity[i][1]=0;

        }

        //      std::cout<<"global Solution of cell"<<i<<"= "
        //           <<this->problem.variables.globalSolution[i]<<std::endl;

  //      std::cout<<"waterDepth = "<<this->problem.variables.waterDepth[i]
    //            <<std::endl;
  //     std::cout<<"velX = "<<this->problem.variables.velocity[i][0]<<std::endl;
        // std::cout<<"velY = "<<this->problem.variables.velocity[i][1]<<std::endl;

    }
    return;
}

}
#endif

// end namespace


