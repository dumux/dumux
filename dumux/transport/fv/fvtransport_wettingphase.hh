// $Id$

#ifndef DUNE_FVTRANSPORT_WETTINGPHASE_HH
#define DUNE_FVTRANSPORT_WETTINGPHASE_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/transport/transport_phasepressure.hh"
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"


/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Markus Wolff, Jochen Fritz
 */

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
     *  of \f$\text{div}\, \boldsymbol{v}_w\f$.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing a CFL condition.
     */
    int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac);

    void initialTransport();

    /*! @brief constructor
     *
     * @param grid a DUNE grid object
     * @param problem an object of class TransportProblem or derived, or a different problem Type
     */
    FVTransport(Grid& grid, Problem& problem)
    :Transport<Grid, Scalar, VC, Problem>(grid, problem)
    {}
};

template<class Grid, class Scalar, class VC, class Problem>
int FVTransport<Grid, Scalar, VC, Problem>::update(const Scalar t, Scalar& dt,
        RepresentationType& updateVec, Scalar& cFLFac = 1)
{
    const GridView& gridView = this->grid_.levelView(this->transProblem.variables.levelTransport);
    // initialize dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec = 0;

    // phase densities
    Scalar densityW = this->transProblem.wettingPhase.density();
    Scalar densityNW = this->transProblem.nonWettingPhase.density();
    Scalar viscosityW = this->transProblem.wettingPhase.viscosity();
    Scalar viscosityNW = this->transProblem.nonWettingPhase.viscosity();
    FieldVector<Scalar,dimWorld> gravity = this->transProblem.gravity();

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
        int globalIdxI = this->transProblem.variables.indexSetTransport.index(*eIt);

        Scalar residualSatW = this->transProblem.soil().Sr_w(globalPos, *eIt, localPos);
        Scalar residualSatNW = this->transProblem.soil().Sr_n(globalPos, *eIt, localPos);
        Scalar porosity = this->transProblem.soil().porosity(globalPos, *eIt,localPos);

        Scalar pressI = this->transProblem.variables.pressure[globalIdxI];
        Scalar pcI = this->transProblem.variables.capillaryPressure[globalIdxI];

        Scalar timestepFactorIn = 0;
        Scalar timestepFactorOutW = 0;
        Scalar timestepFactorOutNW = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator
        isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
            int indexInInside = isIt->indexInInside();

            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition&
            localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(indexInInside,1);

            Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

            // get normal vector scaled with volume
            Dune::FieldVector<Scalar,dimWorld> integrationOuterNormal = isIt->integrationOuterNormal(faceLocal);
            integrationOuterNormal *= Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

            Scalar faceArea = isIt->geometry().volume();

            // get absolute permeability
            FieldMatrix<Scalar,dim,dim> permeabilityI(this->transProblem.soil().K(globalPos, *eIt, localPos));
            // compute directed permeability vector permeabilityI.n
            FieldVector<Scalar,dim> normalPermeabilityI(0);
            permeabilityI.umv(unitOuterNormal, normalPermeabilityI);

            Scalar factor = 0;

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->transProblem.variables.indexSetTransport.index(*neighborPointer);

                // compute flux from one side only
                // this should become easier with the new IntersectionIterator functionality!
                if ( eIt->level()>=neighborPointer->level() )
                {
                    // compute factor in neighbor
                    Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                    const LocalPosition&
                    localPosNeighbor = Dune::ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                    // cell center in global coordinates
                    const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                    // neighbor cell center in global coordinates
                    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosNeighbor - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // get absolute permeability
                    FieldMatrix<Scalar,dim,dim> permeabilityneumann(this->transProblem.soil().K(globalPosNeighbor, *neighborPointer, localPosNeighbor));

                    // compute vectorized permeabilities
                    FieldVector<Scalar,dim> normalPermeabilityneumann(0);
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

                    Scalar pressJ = this->transProblem.variables.pressure[globalIdxJ];
                    Scalar pcJ = this->transProblem.variables.capillaryPressure[globalIdxJ];

                    //calculate velocities
                    Scalar potentialW = (unitOuterNormal * distVec) * (pressI - pressJ) / (dist * dist);
                    Scalar potentialNW = (unitOuterNormal * distVec) * (pressI + pcI - pressJ -pcJ) / (dist * dist);
                     potentialW += densityW * (unitOuterNormal * gravity);
                    potentialNW += densityNW * (unitOuterNormal * gravity);

                    Scalar lambdaW, lambdaNW;

                    if (potentialW >= 0.)
                    {
                        lambdaW = this->transProblem.variables.mobilityWetting[globalIdxI];
                    }
                    else
                    {
                        lambdaW = this->transProblem.variables.mobilityWetting[globalIdxJ];
                    }

                    if (potentialNW >= 0.)
                    {
                        lambdaNW = this->transProblem.variables.mobilityNonWetting[globalIdxI];
                    }
                    else
                    {
                        lambdaNW = this->transProblem.variables.mobilityNonWetting[globalIdxJ];
                    }

                    Scalar velocityW = lambdaW * ((permeability * distVec) * (pressI - pressJ) / (dist * dist) + (permeability * gravity) * densityW);
                     Scalar velocityNW = lambdaNW * ((permeability * distVec) * (pressI + pcI - pressJ - pcJ) / (dist * dist) + (permeability * gravity) * densityW);

                    factor = velocityW * faceArea / (volume*porosity);

                    if (velocityW> 0)
                    {
                        timestepFactorOutW += factor;
                    }
                    if (velocityNW> 0)
                    {
                        timestepFactorOutNW += velocityNW * faceArea / (volume*porosity);
                    }
                    if (velocityW < 0)
                    {
                        Scalar krSum = lambdaW*viscosityW+lambdaNW*viscosityNW;
                        timestepFactorIn -= factor/krSum;
                    }
                    if (velocityNW < 0)
                    {
                        Scalar krSum = lambdaW*viscosityW+lambdaNW*viscosityNW;
                        timestepFactorIn -= velocityNW * faceArea / (volume*porosity*krSum);
                    }
                }
            }

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                //get boundary type
                BoundaryConditions::Flags bcTypeSat = this->transProblem.bctypeSat(globalPosFace, *eIt, localPosFace);
                BoundaryConditions::Flags bcTypePress = this->transProblem.bctypePress(globalPosFace, *eIt, localPosFace);

                // cell center in global coordinates
                GlobalPosition globalPos = eIt->geometry().global(localPos);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                Scalar velocityW = 0;
                Scalar velocityNW = 0;

                Scalar lambdaW = 0, lambdaNW = 0;
                Scalar satBoundEff = 0;
                Scalar satBound = 0;

                if (bcTypeSat == BoundaryConditions::dirichlet)
                {
                    satBound = this->transProblem.dirichletSat(globalPosFace, *eIt, localPosFace);
                    satBoundEff = (satBound-residualSatW)/(1-residualSatW-residualSatNW);
                }
                else
                {
                    satBound = this->transProblem.variables.saturation[globalIdxI];
                    satBoundEff = (satBound -residualSatW)/(1-residualSatW-residualSatNW);
                }

                if (bcTypePress == BoundaryConditions::dirichlet)
                {
                    Scalar pressBound = this->transProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    Scalar pcBound = this->transProblem.materialLaw().pC(satBound, globalPosFace, *eIt, localPosFace);

                    //calculate velocities
                    Scalar potentialW = (unitOuterNormal * distVec) * (pressI - pressBound) / (dist * dist);
                    Scalar potentialNW = (unitOuterNormal * distVec) * (pressI + pcI - pressBound - pcBound) / (dist * dist);

                    potentialW += densityW * (unitOuterNormal * gravity);
                    potentialNW += densityNW * (unitOuterNormal * gravity);

                    if (potentialW >= 0.)
                    {
                        lambdaW = this->transProblem.variables.mobilityWetting[globalIdxI];
                    }
                    else
                    {
                        lambdaW = satBoundEff / viscosityW;
                    }
                    if (potentialNW >= 0.)
                    {
                        lambdaNW = this->transProblem.variables.mobilityNonWetting[globalIdxI];
                    }
                    else
                    {
                        lambdaNW = (1 - satBoundEff) / viscosityNW;
                    }

                    velocityW = lambdaW * ((normalPermeabilityI * distVec) * (pressI - pressBound) / (dist * dist) - (normalPermeabilityI * gravity) * densityW);
                    velocityNW = lambdaNW * ((normalPermeabilityI * distVec) * (pressI + pcI - pressBound - pcBound) / (dist * dist) - (normalPermeabilityI * gravity) * densityW);
                 }
                else
                {
                    Scalar J = this->transProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    if (satBoundEff == 0)
                    {
                        velocityW = 0;
                        velocityNW = J;
                    }
                    else
                    {
                        velocityW=J;
                        velocityNW = 0;
                    }
                }

                factor = velocityW * faceArea / (volume*porosity);

                if (velocityW> 0)
                {
                    timestepFactorOutW += factor;
                }
                if (velocityNW> 0)
                {
                    timestepFactorOutNW += velocityNW * faceArea / (volume*porosity);
                }
                if (velocityW < 0)
                {
                    timestepFactorIn -= factor;
                }
                if (velocityNW < 0)
                {
                    timestepFactorIn -= velocityNW * faceArea / (volume*porosity);
                }

                if (bcTypeSat == BoundaryConditions::neumann)
                {
                    factor = this->transProblem.neumannSat(globalPosFace, *eIt, localPosFace, factor);
                }
            }
            // add to update vector
            updateVec[globalIdxI] -= factor;
        }
        // end all intersections
        // compute dt restriction
        Scalar volumeCorrectionFactorIn = (1-residualSatW - residualSatNW);
        Scalar volumeCorrectionFactorOutW = (this->transProblem.variables.saturation[globalIdxI]-residualSatW);
        Scalar volumeCorrectionFactorOutNW = (1-this->transProblem.variables.saturation[globalIdxI] - residualSatNW);

        //make sure correction is in the right range. If not: force dt to be not min-dt!
        if (volumeCorrectionFactorOutW <= 0)
        {
            volumeCorrectionFactorOutW = 1e100;
        }
        if (volumeCorrectionFactorOutNW <= 0)
        {
            volumeCorrectionFactorOutNW = 1e100;
        }

        //make sure correction is in the right range. If not: force dt to be not min-dt!
        if (timestepFactorIn <= 0)
        {
            timestepFactorIn = 1e-100;
        }
        if (timestepFactorOutW <= 0)
        {
            timestepFactorOutW = 1e-100;
        }
        if (timestepFactorOutNW <= 0)
        {
            timestepFactorOutNW = 1e-100;
        }

        timestepFactorIn = volumeCorrectionFactorIn/timestepFactorIn;
        timestepFactorOutW = volumeCorrectionFactorOutW/timestepFactorOutW;
        timestepFactorOutNW = volumeCorrectionFactorOutNW/timestepFactorOutNW;

        Scalar timestepFactor = std::min(timestepFactorOutW,timestepFactorOutNW);
        timestepFactor = std::min(timestepFactor, timestepFactorIn);

        dt = std::min(dt, timestepFactor);

    } // end grid traversal

    return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
void FVTransport<Grid, Scalar, VC, Problem>::initialTransport()
{
    const GridView& gridView = this->grid_.levelView(this->transProblem.variables.levelTransport);
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
        this->transProblem.variables.saturation[this->transProblem.variables.indexSetTransport.index(*eIt)] = this->transProblem.initSat(globalPos, *eIt, localPos);
    }
    return;
}

}
#endif
