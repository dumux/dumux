// $Id$

#ifndef DUNE_FVTRANSPORT_HH
#define DUNE_FVTRANSPORT_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/transport/transport.hh"
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transport
//! The finite volume model for the solution of the transport equation
template<class Grid, class Scalar, class VC, class Problem = TransportProblem<
        Grid, Scalar, VC> >
class FVTransport: public Transport<Grid, Scalar, VC, Problem>
{
    enum
    {
        dim = Grid::dimension
    };
    enum
    {
        dimWorld = Grid::dimensionworld
    };
typedef    typename VC::ScalarVectorType PressType;
    typedef typename VC::VelType VelType;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef BlockVector< Dune::FieldVector<Scalar,dim> > SlopeType;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef typename VC::ScalarVectorType RepresentationType;
    /*!
     *  \param t time
     *  \param dt time step size to estimate
     *  \param update vector to be filled with the update
     *
     *  This method calculates the update vector, i.e., the FV discretization
     *  of \f$\text{div}\, (f_\text{w} \boldsymbol{v}_t)\f$. The total velocity
     *  \f$\boldsymbol{v}_t)\f$ has to be given by the block vector \a velocity,
     *  containing values at each cell face. The fractional flow function \f$f_\text{w}\f$
     *  is realized by the numerical flux function \a numericalFlux such that the \f$i\f$-th entry
     *  of the vector update is obtained by calculating
     *  \f[ \sum_{j \in \mathcal{N}(i)} v_{ij}g(S_i, S_j) - v_{ji}g(S_j, S_i), \f]
     *  where \f$\mathcal{N}(i)\f$ denotes the set of neighbors of the cell \f$i\f$ and
     *  \f$g\f$ stands for \a numericalFlux. The normal velocities \f$v_{ij}\f$ and \f$v_{ji}\f$
     *  are given by
     *  \f{align*} v_{ij} = \int_{\Gamma_{ij}} \max(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \qquad
     *  v_{ji} = \int_{\Gamma_{ij}} \min(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \f}
     *  where \f$\boldsymbol{n}_{ij}\f$ denotes the unit normal vector from cell \f$i\f$ towards cell \f$j\f$.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing the usual CFL condition.
     */
    int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac);

    void initialTransport();

    /*! @brief constructor
     *
     * @param g a DUNE grid object
     * @param prob an object of class TransportProblem or derived
     * @param lev the grid level on which the Transport equation is to be solved.
     * @param diffPart an object of class DiffusivePart or derived. This determines the diffusive flux incorporated in the transport.
     * @param rec flag to switch on linear reconstruction (second order TVD)
     * @param amax alphamax parameter for slope limiter in TVD
     * @param numFl an object of class Numerical Flux or derived
     */
    FVTransport(Grid& grid, Problem& problem,
            DiffusivePart<Grid,Scalar>& diffPart = *(new DiffusivePart<Grid, Scalar>), bool rec = false,
            Scalar aMax = 0.8, const NumericalFlux<Scalar>& numFl = *(new Upwind<Scalar>)) :
    Transport<Grid, Scalar, VC, Problem>(grid, problem),
    reconstruct_(rec),numFlux_(numFl), diffusivePart_(diffPart), alphaMax_(aMax)
    {}

private:

    void calculateSlopes(SlopeType& slope, Scalar t, Scalar& cFLFactor);

    bool reconstruct_;
    const NumericalFlux<Scalar>& numFlux_;
    const DiffusivePart<Grid, Scalar>& diffusivePart_;
    Scalar alphaMax_;
};

template<class Grid, class Scalar, class VC, class Problem>
int FVTransport<Grid, Scalar, VC, Problem>::update(const Scalar t, Scalar& dt,
        RepresentationType& updateVec, Scalar& cFLFac = 1)
{
    const GridView& gridView = this->grid_.levelView(this->level());
    // initialize dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec = 0;

    SlopeType slope(this->transProblem.variables.transMapper.size());
    calculateSlopes(slope, t,cFLFac);

    // compute update vector
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        //
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)
        *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        // cell index
        int globalIdxI = this->transProblem.variables.transMapper.map(*eIt);

        // for time step calculation
        Scalar sumFactor = 0;
        Scalar sumFactor2 = 0;
        Scalar sumDiff = 0;
        Scalar sumDiff2 = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator
        isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
            int numberInSelf = isIt->indexInInside();

            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition&
            localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(numberInSelf,1);

            // get normal vector scaled with volume
            Dune::FieldVector<Scalar,dimWorld> integrationOuterNormal = isIt->integrationOuterNormal(faceLocal);
            integrationOuterNormal *= Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

            // compute factor occuring in flux formula
            Scalar velocityIJ = std::max(this->transProblem.variables.vTotal(*eIt, numberInSelf)*integrationOuterNormal/(volume), 0.0);

            Scalar factor = 0, diffFactor = 0, totFactor = 0;

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->transProblem.variables.transMapper.map(*neighborPointer);

                // compute flux from one side only
                // this should become easier with the new IntersectionIterator functionality!
                if ( eIt->level()>=neighborPointer->level() )
                {
                    // compute factor in neighbor
                    Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                    const LocalPosition&
                    localPosNeighbor = Dune::ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                    Scalar velocityJI = std::max(-(this->transProblem.variables.vTotal(*eIt, numberInSelf)*integrationOuterNormal/volume), 0.0);

                    // cell center in global coordinates
                    const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                    // neighbor cell center in global coordinates
                    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosNeighbor;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // get saturation value at cell center
                    Scalar satI = this->transProblem.variables.saturation[globalIdxI];

                    // get saturation value at neighbor cell center
                    Scalar satJ = this->transProblem.variables.saturation[globalIdxJ];

                    // calculate the saturation gradient
                    Dune::FieldVector<Scalar,dim> satGradient = distVec;
                    satGradient *= (satJ - satI)/(dist*dist);

                    // the arithmetic average
                    Scalar satAvg = 0.5*(satI + satJ);

                    // get the diffusive part
                    Scalar diffPart = diffusivePart_(*eIt, numberInSelf, satAvg, satGradient, t, satI, satJ)*integrationOuterNormal;

                    // CAREFUL: works only for axisymmetric grids
                    if (reconstruct_)
                    {
                        for (int k = 0; k < dim; k++)
                        if (fabs(distVec[k])> 0.5*dist)
                        {
                            satI -= fabs(distVec[k])/distVec[k]*0.5*dist*slope[globalIdxI][k];
                            satJ += fabs(distVec[k])/distVec[k]*0.5*dist*slope[globalIdxJ][k];
                        }
                    }

                    Scalar fI = this->transProblem.materialLaw().fractionalW(satI, globalPos, *eIt,localPos);
                    Scalar fJ = this->transProblem.materialLaw().fractionalW(satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor);

                    diffFactor = diffPart / volume;
                    factor = diffFactor
                    + velocityJI*numFlux_(satJ, satI, fJ, fI)
                    - velocityIJ*numFlux_(satI, satJ, fI, fJ);
                    factor/= this->transProblem.soil().porosity(globalPos, *eIt,localPos);

                    std::vector<Scalar> relativePermeabilitiesI = this->transProblem.materialLaw().kr(satI, globalPos, *eIt,localPos);
                    std::vector<Scalar> relativePermeabilitiesJ = this->transProblem.materialLaw().kr(satJ, globalPosNeighbor, *neighborPointer,localPosNeighbor);

                    Scalar krSumI = relativePermeabilitiesI[0] + relativePermeabilitiesI[1];
                    Scalar krSumJ = relativePermeabilitiesJ[0] + relativePermeabilitiesJ[1];

                    totFactor = velocityJI/krSumJ - velocityIJ/krSumI;
                }
            }

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                //get boundary type
                BoundaryConditions::Flags bctype = this->transProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    // get saturation value at cell center
                    Scalar satI = this->transProblem.variables.saturation[globalIdxI];

                    Scalar velocityJI = std::max(-(this->transProblem.variables.vTotal(*eIt, numberInSelf)*integrationOuterNormal/volume), 0.0);

                    Scalar satBound = this->transProblem.dirichletSat(globalPosFace, *eIt, localPosFace);

                    // cell center in global coordinates
                    GlobalPosition globalPos = eIt->geometry().global(localPos);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // calculate the saturation gradient
                    Dune::FieldVector<Scalar,dim> satGradient = distVec;
                    satGradient *= (satBound - satI)/(dist*dist);

                    // the arithmetic average
                    Scalar satAvg = 0.5*(satI + satBound);

                    // get the diffusive part
                    Scalar diffPart = diffusivePart_(*eIt, numberInSelf, satAvg, satGradient, t, satI, satBound)*integrationOuterNormal;

                    // CAREFUL: works only for axisymmetric grids
                    if (reconstruct_)
                    {
                        for (int k = 0; k < dim; k++)
                        if (fabs(distVec[k])> 0.5*dist)
                        {
                            //TODO remove DEBUG--->
                            //Scalar gabagabahey = slope[globalIdxI][k];
                            //<---DEBUG
                            satI -= fabs(distVec[k])/distVec[k]*dist*slope[globalIdxI][k];
                        }
                    }

                    Scalar fI = this->transProblem.materialLaw().fractionalW(satI, globalPos, *eIt,localPos);
                    Scalar fBound = this->transProblem.materialLaw().fractionalW(satBound, globalPosFace, *eIt, localPosFace);
                    diffFactor = diffPart / volume;
                    factor = diffFactor
                    + velocityJI*numFlux_(satBound, satI, fBound, fI)
                    - velocityIJ*numFlux_(satI, satBound, fI, fBound);
                    factor/= this->transProblem.soil().porosity(globalPos, *eIt,localPos);

                    std::vector<Scalar> relativePermeabilitiesI = this->transProblem.materialLaw().kr(satI, globalPos, *eIt,localPos);
                    std::vector<Scalar> relativePermeabilitiesBound = this->transProblem.materialLaw().kr(satBound, globalPosFace, *eIt,localPosFace);

                    Scalar krSumI = relativePermeabilitiesI[0] + relativePermeabilitiesI[1];
                    Scalar krSumBound = relativePermeabilitiesBound[0] + relativePermeabilitiesBound[1];

                    totFactor = velocityJI/krSumBound - velocityIJ/krSumI;
                }
                else
                {
                    Scalar satI = this->transProblem.variables.saturation[globalIdxI];
                    Scalar velocityJI = std::max(-(this->transProblem.variables.vTotal(*eIt, numberInSelf)*integrationOuterNormal/volume), 0.0);
                    Scalar fI = this->transProblem.materialLaw().fractionalW(satI, globalPos, *eIt,localPos);
                    Scalar helpFactor = (velocityJI-velocityIJ)*fI;
                    factor = this->transProblem.neumannSat(globalPosFace, *eIt, localPosFace, helpFactor);
                    factor/= this->transProblem.soil().porosity(globalPos, *eIt,localPos);

                    totFactor = 0;
                    if (factor)
                    {
                        //                        std::cout<<"factor = "<<factor<<std::endl;
                        //                        std::cout<<"updateVec = "<<updateVec[globalIdxI]<<std::endl;
                        totFactor = velocityJI-velocityIJ;
                    }
                }
            }
            // add to update vector
            updateVec[globalIdxI] += factor;

            // for time step calculation
            if (totFactor>=0)
            sumFactor += totFactor;
            else
            sumFactor2 += (-totFactor);
            if (diffFactor>=0)
            sumDiff += diffFactor;
            else
            sumDiff2 += (-diffFactor);
        }
        Scalar volumeCorrectionFactor = (1-this->transProblem.soil().Sr_w(globalPos, *eIt,localPos)-this->transProblem.soil().Sr_n(globalPos, *eIt,localPos))*this->transProblem.soil().porosity(globalPos, *eIt,localPos);

        // end all intersections
        // compute dt restriction
        sumFactor = std::min(volumeCorrectionFactor/sumFactor, volumeCorrectionFactor/sumFactor2);
        if (sumDiff> 0 && sumDiff2> 0)
        {
            sumDiff = std::min(volumeCorrectionFactor/sumDiff, volumeCorrectionFactor/sumDiff2);
            sumFactor = std::min(sumFactor, 0.1*sumDiff);
        }
        dt = std::min(dt, sumFactor);

    } // end grid traversal

    return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
void FVTransport<Grid, Scalar, VC, Problem>::initialTransport()
{
    const GridView& gridView = this->grid_.levelView(this->level());
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // initialize cell concentration
        this->transProblem.variables.saturation[this->transProblem.variables.transMapper.map(*eIt)] = this->transProblem.initSat(globalPos, *eIt, localPos);
    }
    return;
}

template<class Grid, class Scalar, class VC, class Problem>
void FVTransport<Grid, Scalar, VC, Problem>::calculateSlopes(
        SlopeType& slope, Scalar t, Scalar& cFLFactor)
{
    const GridView& gridView = this->grid_.levelView(this->level());
    Scalar stabilityFactor = 1.0 - cFLFactor*sqrt(cFLFactor);

    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt!=eItEnd; ++eIt)
    {
        // get some cell properties
        Dune::GeometryType gt = eIt->geometry().type();
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);
        GlobalPosition globalPos = eIt->geometry().global(localPos);
        int globalIdxI = this->transProblem.variables.transMapper.map(*eIt);

        // vector containing the distances to the neighboring cells
        Dune::FieldVector<Scalar, 2*dim> dist;

        // vector containing the saturations of the neighboring cells
        Dune::FieldVector<Scalar, 2*dim> saturation;

        // location[k], k = 0,...,2dim-1, contains the local index w.r.t. IntersectionIterator,
        // i.e. the numberInSelf, which isIt east, west, north, south, top, bottom to the cell center
        Dune::FieldVector<int, 2*dim> location;

        // run through all intersections with neighbors and boundary
        IntersectionIterator
        isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
            int numberInSelf = isIt->indexInInside();

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->transProblem.variables.transMapper.map(*neighborPointer);

                // get saturation value
                saturation[numberInSelf] = this->transProblem.variables.saturation[globalIdxJ];

                // compute factor in neighbor
                Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition
                &localPosNeighbor = Dune::ReferenceElements<Scalar,dim>::general(neighborGT).position(0, 0);

                // neighboring cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosNeighbor;

                // compute distance between cell centers
                dist[numberInSelf] = distVec.two_norm();

                // CAREFUL: works only for axiparallel grids
                for (int k = 0; k < dim; k++)
                if (globalPosNeighbor[k] - globalPos[k]> 0.5*dist[numberInSelf])
                {
                    location[2*k] = numberInSelf;
                }
                else if (globalPosNeighbor[k] - globalPos[k]< -0.5*dist[numberInSelf])
                {
                    location[2*k + 1] = numberInSelf;
                }
            }

            // handle boundary face
            if (isIt->boundary())
            {
                // get geometry type of face
                Dune::GeometryType faceGT = isIt->geometryInInside().type();

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>&
                faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                // center of face in global coordinates
                const GlobalPosition&
                globalPosFace = isIt->geometry().global(faceLocal);

                // get saturation value
                saturation[numberInSelf] = this->transProblem.variables.saturation[globalIdxI];

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;

                // compute distance
                dist[numberInSelf] = distVec.two_norm();

                // CAREFUL: works only for axiparallel grids
                for (int k = 0; k < dim; k++)
                if (globalPosFace[k] - globalPos[k]> 0.5*dist[numberInSelf])
                {
                    location[2*k] = numberInSelf;
                }
                else if (globalPosFace[k] - globalPos[k] < -0.5*dist[numberInSelf])
                {
                    location[2*k + 1] = numberInSelf;
                }
            }
        } // end all intersections

        for (int k = 0; k < dim; k++)
        {
            Scalar slopeIK = (saturation[location[2*k]]
                    - saturation[location[2*k + 1]])/(dist[location[2*k]]
                    + dist[location[2*k + 1]]);

            Scalar alphaIK = 1.0;
            if (fabs(slopeIK)> 1e-8*dist[location[2*k]])
            {
                Scalar satI = this->transProblem.variables.saturation[globalIdxI];
                Scalar alphaRight = stabilityFactor*fabs(2.0
                        /(dist[location[2*k]]*slopeIK)
                        *(saturation[location[2*k]] - satI));
                Scalar alphaLeft = stabilityFactor*fabs(2.0
                        /(dist[location[2*k + 1]]*slopeIK)*(satI
                                - saturation[location[2*k + 1]]));
                alphaIK = std::min(std::min(std::max(alphaRight, 0.0),
                                std::max(alphaLeft, 0.0)), alphaMax_);
            }

            slope[globalIdxI][k] = alphaIK*slopeIK;

        }
    } // end grid traversal

    return;
}

}
#endif
