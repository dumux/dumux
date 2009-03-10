// $Id:$

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
//! The finite volume model for the solution of the transport equation
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

    int update(const Scalar t, Scalar& dt, SolutionType& updateVec,
               Scalar& cFLFac);

    void initialize();

    FVShallowWater(Grid& grid, ShallowProblemBase<Grid, Scalar, VC>& problem,
                   NumericalFlux<Grid,Scalar>& numFl = *(new FirstOrderUpwind<Grid,Scalar>)) :
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
                                                                VC>::update(const Scalar t, Scalar& dt, SolutionType& updateVec,
                                                                            Scalar& cFLFac = 1)
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
        //Scalar dist;
        Scalar wDepthI;
        Scalar wDepthJ;
        Scalar wDepthFace;
        //Scalar wDepthFaceI;
        //Scalar wDepthFaceJ;
        SlopeType velI;
        SlopeType velJ;
        SlopeType bottomSlopeTerm(0); //slope term resulting from divergence form
        SlopeType bottomSlope(0); //real bottomslope
        SystemType bottomSlopeVector(0);
        SystemType flux(0);
        SlopeType divergenceTerm(0);
        SystemType summedFluxes(0);
        Scalar gravity(0);
        SlopeType velFace(0);

        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition &localPos =
            Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)
            *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        //get bottomslopevector for entity
        bottomSlope =this->problem.surface.calcBottomSlopes(globalPos, *eIt,
                                                            localPos);

        // get waterdepth at cell center
        wDepthI = this->problem.variables.wDepth[globalIdxI];

        // get velocity at cell center
        velI = this->problem.variables.velocity[globalIdxI];

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd =gridView.template iend(*eIt);
        for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
        {

            // local number of facet
            // int numberInSelf = isIt->numberInSelf();

            // get geometry type of face
            Dune::GeometryType gtf = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>& faceLocalPos =
                Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0, 0);

            // center of face inside volume reference element
            const LocalPosition& faceLocalDim =
                Dune::ReferenceElements<Scalar,dim>::general(gtf).position(isIt->numberInSelf(), 1);

            //get normal vector of face
            Dune::FieldVector<Scalar,dimWorld> nVec =
                isIt->outerNormal(faceLocalPos);

            // get normal vector scaled with volume
            Dune::FieldVector<Scalar,dimWorld> nVecScaled =
                isIt->integrationOuterNormal(faceLocalPos);
            nVecScaled*=Dune::ReferenceElements<Scalar,dim-1>::general(gtf).volume();

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
                //Scalar dist = distVec.two_norm();

                // get waterdepth at neighbor cell center
                wDepthJ = this->problem.variables.wDepth[globalIdxJ];

                // get velocity at neighbor cell center
                velJ = this->problem.variables.velocity[globalIdxJ];

                //fluxVector for a given scheme (Upwindor different will be catched from numericalflux)

                flux = numFlux(velI, velJ, wDepthI, wDepthJ, nVecScaled, nVec);

                // std::cout<<"cell "<<globalIdxI<<"flux_interior of face "
                //<<numberInSelf<<"=" <<flux<<std::endl;

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
                //Scalar dist = distVec.two_norm();

                //get boundary type
                BoundaryConditions::Flags bctype = this->problem.bctype(
                                                                        faceGlobalPos, *eIt, faceLocalDim);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    // get waterdepth at boundary
                    wDepthFace = this->problem.dirichlet(faceGlobalPos, *eIt,
                                                         faceLocalDim);

                    Scalar conti = velFace * wDepthFace;
                    Scalar xMomentum = velFace[0] * wDepthFace*velFace;
                    xMomentum += gravity*wDepthFace*0.5;

                    Scalar yMomentum = velFace[1] * wDepthFace*velFace;
                    xMomentum += gravity*wDepthFace*0.5;

                    flux[0] = conti;
                    flux[1] = xMomentum;
                    flux[2] = yMomentum;

                    //std::cout<<"flux_dirichlet of face "<<numberInSelf<<"="
                    //<<flux<<std::endl;


                }
                if (bctype == BoundaryConditions::neumann) //no flow condition at side walls

                {

                    flux = this->problem.neumann(faceGlobalPos, *eIt,
                                                 faceLocalDim);

                    //std::cout<<"flux_neumann of face "<<numberInSelf<<"="<<flux
                    //<<std::endl;
                }
            }
            summedFluxes += flux;

        }

        bottomSlopeVector[0]=0;
        bottomSlopeVector[1]=9.81;
        bottomSlopeVector[1]*=bottomSlope[0];
        bottomSlopeVector[1]*=(wDepthI);
        bottomSlopeVector[2]=9.81;
        bottomSlopeVector[2]*=bottomSlope[1];
        bottomSlopeVector[2]*=(wDepthI);

        //std::cout<<"BottomSlope Vector of cell="<<globalIdxI<<"= "
        //<<bottomSlopeVector<<std::endl;

        //Quelle einbinden und in Systemvektor konvertieren
        Scalar sourceTerm = this->problem.setSource(globalPos, *eIt, localPos);

        SystemType sourceTermVector(0);
        sourceTermVector[0]=sourceTerm;
        sourceTermVector[1]=0;
        sourceTermVector[2]=0;
        sourceTermVector *= volume;

        //std::cout<<"Source Term Vector of cell="<<globalIdxI<<"= "
        //<<sourceTermVector<<std::endl;

        updateVec[globalIdxI]=summedFluxes;
        updateVec[globalIdxI]-= sourceTermVector;
        updateVec[globalIdxI]-=bottomSlopeVector;
        updateVec[globalIdxI]/= (-volume);

        //std::cout<<"update for cell"<<globalIdxI<<" = "<<updateVec[globalIdxI]
        //<<std::endl;

        //add cfl-criterium
        //calculate timestep for every cell and take the minimum

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

        Scalar initialWaterDepth=this->problem.setInitWDepth(globalPos, *eIt,
                                                             localPos);
        SlopeType initialVelocity=this->problem.setInitVel(globalPos, *eIt,
                                                           localPos);

        // initialize cell values
        this->problem.variables.wDepth[globalIdx] = initialWaterDepth;
        this->problem.variables.velocity[globalIdx] = initialVelocity;
        this->problem.variables.globalSolution[globalIdx][0]
            = initialWaterDepth;
        this->problem.variables.globalSolution[globalIdx][1]
            = initialWaterDepth*initialVelocity[0];
        this->problem.variables.globalSolution[globalIdx][1]
            = initialWaterDepth*initialVelocity[1];
    }
    return;
}

}
#endif

// end namespace


