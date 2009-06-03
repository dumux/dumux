// $Id$

#ifndef DUNE_FVSATURATIONWETTING2P_HH
#define DUNE_FVSATURATIONWETTING2P_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/transport/transport.hh"
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
template<class GridView, class Scalar, class VC,
        class Problem = TransportProblem<GridView, Scalar, VC> >
class FVSaturationWetting2P: public Transport<GridView, Scalar, VC, Problem>
{
    enum
    {
        dim = GridView::dimension
    };
    enum
    {
        dimWorld = GridView::dimensionworld
    };
    enum
    {
        vw = 0, vt = 1
    };
typedef    typename VC::ScalarVectorType PressType;
    typedef typename VC::VelType VelType;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
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
    int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac, bool impes);

    void initialTransport();

    Scalar evaluateTimeStepWettingFlux(Scalar timestepFactorIn, Scalar timestepFactorOutW ,Scalar& residualSatW, Scalar& residualSatNW, int globalIdxI);

    Scalar evaluateTimeStepTotalFlux(Scalar timestepFactorIn,Scalar timestepFactorOut, Scalar diffFactorIn, Scalar diffFactorOut, Scalar& residualSatW, Scalar& residualSatNW);

    void updateMaterialLaws();

    virtual void vtkout(const char* name, int k) const
    {
        this->transProblem.variables().vtkout(name, k);
        return;
    }

    /*! @brief constructor
     *
     * @param grid a DUNE grid object
     * @param problem an object of class TransportProblem or derived, or a different problem Type
     */
    FVSaturationWetting2P(GridView& gridView, Problem& problem, std::string velocityType, DiffusivePart<GridView,Scalar>& diffPart = *(new DiffusivePart<GridView, Scalar>))
    :Transport<GridView, Scalar, VC, Problem>(gridView, problem), diffusivePart_(diffPart),
    velocityType_((velocityType == "vw") ? 0 : ((velocityType == "vt") ? 1 : 999))
    {
        if (velocityType_ == 999)
        {
            DUNE_THROW(NotImplemented, "Velocity type not supported!");
        }
    }
private:
    const DiffusivePart<GridView, Scalar>& diffusivePart_;
    const int velocityType_;
};

template<class GridView, class Scalar, class VC, class Problem>
int FVSaturationWetting2P<GridView, Scalar, VC, Problem>::update(const Scalar t, Scalar& dt,
        RepresentationType& updateVec, Scalar& cFLFac = 1, bool impes = false)
{
    if (!impes)
    {
        updateMaterialLaws();
    }

    // initialize dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec = 0;

    // some phase properties
    Scalar viscosityW = this->transProblem.wettingPhase().viscosity();
    Scalar viscosityNW = this->transProblem.nonWettingPhase().viscosity();
    Scalar viscosityRatio = 1 - fabs(viscosityW-viscosityNW)/(viscosityW+viscosityNW);
    FieldVector<Scalar,dimWorld> gravity = this->transProblem.gravity();

    // compute update vector
    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
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
        int globalIdxI = this->transProblem.variables().indexTransport(*eIt);

        Scalar residualSatW = this->transProblem.soil().Sr_w(globalPos, *eIt, localPos);
        Scalar residualSatNW = this->transProblem.soil().Sr_n(globalPos, *eIt, localPos);
        Scalar porosity = this->transProblem.soil().porosity(globalPos, *eIt,localPos);

        Scalar timestepFactorIn = 0;
        Scalar timestepFactorOut = 0;
        Scalar timestepFactorOutW = 0;
        Scalar diffFactorIn = 0;
        Scalar diffFactorOut = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator
        isItEnd = this->gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = this->gridView.template ibegin(*eIt); isIt
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

            Dune::FieldVector<Scalar,dimWorld> unitOuterNormalAbs(unitOuterNormal);
            for (int i=0;i<dim;i++)
            {
                unitOuterNormalAbs[i] *= unitOuterNormal[i];
            }

            Scalar faceArea = isIt->geometry().volume();

            Scalar factor = 0;

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->transProblem.variables().indexTransport(*neighborPointer);

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

                Scalar potentialW = this->transProblem.variables().potentialWetting()[globalIdxI][indexInInside];
                Scalar potentialNW = this->transProblem.variables().potentialNonWetting()[globalIdxI][indexInInside];

                factor = (this->transProblem.variables().velocity()[globalIdxI][indexInInside] * unitOuterNormal) * faceArea / (volume*porosity);

                if (velocityType_ == vt)
                {
                    Scalar lambdaW, lambdaNW;

                    if (potentialW >= 0.)
                    {
                        lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                    }
                    else
                    {
                        lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxJ];
                    }

                    if (potentialNW >= 0.)
                    {
                        lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                    }
                    else
                    {
                        lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxJ];
                    }

                    //for time step criterion
                    Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                    if (factor >= 0)
                    {
                        timestepFactorOut += factor/(krSum*viscosityRatio);
                    }
                    if (factor < 0)
                    {
                        timestepFactorIn -= factor/(krSum*viscosityRatio);
                    }

                    Scalar satI = this->transProblem.variables().saturation()[globalIdxI];
                    Scalar satJ = this->transProblem.variables().saturation()[globalIdxJ];

                    Scalar pcI = this->transProblem.variables().capillaryPressure()[globalIdxI];
                    Scalar pcJ = this->transProblem.variables().capillaryPressure()[globalIdxJ];

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // calculate the capillary gradient
                    Dune::FieldVector<Scalar,dimWorld> pcGradient = distVec;
                    pcGradient *= (pcJ - pcI)/(dist*dist);

                    // get the diffusive part
                    Scalar diffPart = diffusivePart_(*eIt, indexInInside, satI, satJ, pcGradient)*unitOuterNormal * faceArea / (volume*porosity);

                    if (diffPart >= 0)
                    {
                        diffFactorOut += diffPart/(krSum*viscosityRatio);
                    }
                    if (diffPart < 0)
                    {
                        diffFactorIn -= diffPart/(krSum*viscosityRatio);
                    }

                    //vt*fw
                    factor *= lambdaW/(lambdaW + lambdaNW);
                    factor += diffPart;
                }

                //for time step criterion
                if (velocityType_ == vw)
                {
                    if (potentialW >= 0)
                    {
                        timestepFactorOut+= factor;

                        if (potentialNW < 0)
                        {
                            Scalar lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                            Scalar lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxJ];
                            Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                            //                            std::cout<<(factor/(lambdaNW*potentialNW)*lambdaW*potentialW)<<std::endl;
                            timestepFactorIn -= factor/(lambdaW*potentialW*krSum*viscosityRatio)*lambdaNW*potentialNW;
                        }

                    }
                    if (potentialW < 0)
                    {
                        Scalar lambdaNW = 0;
                        Scalar lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxJ];
                        Scalar krSum = 1;
                        if (potentialNW >= 0)
                        {
                            lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                            krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                        }
                        if (potentialNW < 0)
                        {
                            lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxJ];
                            krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                            timestepFactorIn -= factor/(lambdaW*potentialW*krSum*viscosityRatio)*lambdaNW*potentialNW;
                        }

                        timestepFactorIn -= factor/(krSum*viscosityRatio);
                    }
                    if (isnan(timestepFactorIn) || isinf(timestepFactorIn))
                    {
                        timestepFactorIn = 1e-100;
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

                if (bcTypeSat == BoundaryConditions::dirichlet)
                {
                    // cell center in global coordinates
                    GlobalPosition globalPos = eIt->geometry().global(localPos);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                    Scalar satBound = this->transProblem.dirichletSat(globalPosFace, *eIt, localPosFace);

                    Scalar potentialW = this->transProblem.variables().potentialWetting()[globalIdxI][indexInInside];
                    Scalar potentialNW = this->transProblem.variables().potentialNonWetting()[globalIdxI][indexInInside];

                    factor = (this->transProblem.variables().velocity()[globalIdxI][indexInInside] * unitOuterNormal) * faceArea / (volume*porosity);

                    if (velocityType_ == vt)
                    {
                        Scalar lambdaW, lambdaNW;

                        if (potentialW >= 0.)
                        {
                            lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                        }
                        else
                        {
                            lambdaW = this->transProblem.materialLaw().mobW(satBound,globalPosFace, *eIt, localPosFace);
                        }

                        if (potentialNW >= 0.)
                        {
                            lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                        }
                        else
                        {
                            lambdaNW = this->transProblem.materialLaw().mobN((1-satBound),globalPosFace, *eIt, localPosFace);
                        }

                        //for time step criterion
                        Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                        if (factor >= 0)
                        {
                            timestepFactorOut += factor/(krSum*viscosityRatio);
                        }
                        if (factor < 0)
                        {
                            timestepFactorIn -= factor/(krSum*viscosityRatio);
                        }

                        Scalar satI = this->transProblem.variables().saturation()[globalIdxI];

                        Scalar pcI = this->transProblem.variables().capillaryPressure()[globalIdxI];
                        Scalar pcBound = this->transProblem.materialLaw().pC(satBound,globalPosFace, *eIt, localPosFace);

                        // compute distance between cell centers
                        Scalar dist = distVec.two_norm();

                        // calculate the capillary gradient
                        Dune::FieldVector<Scalar,dimWorld> pcGradient = distVec;
                        pcGradient *= (pcBound - pcI)/(dist*dist);

                        // get the diffusive part
                        Scalar diffPart = diffusivePart_(*eIt, indexInInside, satI, satBound, pcGradient)* unitOuterNormal * faceArea / (volume*porosity);

                        if (diffPart >= 0)
                        {
                            diffFactorOut += diffPart/(krSum*viscosityRatio);
                        }
                        if (diffPart < 0)
                        {
                            diffFactorIn -= diffPart/(krSum*viscosityRatio);
                        }

                        //vt*fw
                        factor *= lambdaW/(lambdaW + lambdaNW);
                        factor += diffPart;
                    }

                    //for time step criterion
                    if (velocityType_ == vw)
                    {
                        if (potentialW >= 0)
                        {
                            timestepFactorOut+= factor;

                            if (potentialNW < 0)
                            {
                                Scalar lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                                Scalar lambdaNW = this->transProblem.materialLaw().mobN(1-satBound,globalPosFace, *eIt, localPosFace);
                                Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                                timestepFactorIn -= factor/(lambdaW*potentialW*krSum*viscosityRatio)*lambdaNW*potentialNW;
                            }

                        }
                        if (potentialW < 0)
                        {
                            Scalar lambdaNW = 0;
                            Scalar lambdaW = this->transProblem.materialLaw().mobW(satBound,globalPosFace, *eIt, localPosFace);
                            Scalar krSum = 1;
                            if (potentialNW >= 0)
                            {
                                lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                                krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                            }
                            if (potentialNW < 0)
                            {
                                lambdaNW = this->transProblem.materialLaw().mobW(1-satBound, globalPosFace, *eIt, localPosFace);
                                krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                                timestepFactorIn -= factor/(lambdaW*potentialW*krSum*viscosityRatio)*lambdaNW*potentialNW;
                            }

                            timestepFactorIn -= factor/(krSum*viscosityRatio);
                        }
                        if (isnan(timestepFactorIn) || isinf(timestepFactorIn))
                        {
                            timestepFactorIn = 1e-100;
                        }
                    }
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
        if (velocityType_ == vw)
        {
            dt = std::min(dt, evaluateTimeStepWettingFlux(timestepFactorIn, timestepFactorOutW , residualSatW, residualSatNW, globalIdxI));
        }
        if (velocityType_ == vt)
        {
            dt = std::min(dt, evaluateTimeStepTotalFlux(timestepFactorIn, timestepFactorOut, diffFactorIn, diffFactorOut, residualSatW, residualSatNW));
        }
    } // end grid traversal

    return 0;
}

template<class GridView, class Scalar, class VC, class Problem>
void FVSaturationWetting2P<GridView, Scalar, VC, Problem>::initialTransport()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // initialize cell concentration
        this->transProblem.variables().saturation()[this->transProblem.variables().indexTransport(*eIt)] = this->transProblem.initSat(globalPos, *eIt, localPos);
    }

    return;
}

template<class GridView, class Scalar, class VC, class Problem>
Scalar FVSaturationWetting2P<GridView, Scalar, VC, Problem>::evaluateTimeStepTotalFlux(Scalar timestepFactorIn,Scalar timestepFactorOut, Scalar diffFactorIn, Scalar diffFactorOut, Scalar& residualSatW, Scalar& residualSatNW)
{
    Scalar volumeCorrectionFactor = (1 - residualSatW -residualSatNW);

    if (timestepFactorIn <= 0)
    {
        timestepFactorIn = 1e-100;
    }
    if (timestepFactorOut <= 0)
    {
        timestepFactorOut = 1e-100;
    }

    Scalar sumFactor = std::min(volumeCorrectionFactor/timestepFactorIn, volumeCorrectionFactor/timestepFactorOut);

    if (diffFactorIn <= 0)
    {
        diffFactorIn = 1e-100;
    }
    if (diffFactorOut <= 0)
    {
        diffFactorOut = 1e-100;
    }

    Scalar minDiff = std::min(volumeCorrectionFactor/diffFactorIn,volumeCorrectionFactor/diffFactorOut);

    sumFactor = std::min(sumFactor, 0.1*minDiff);

    return sumFactor;
}
template<class GridView, class Scalar, class VC, class Problem>
Scalar FVSaturationWetting2P<GridView, Scalar, VC, Problem>::evaluateTimeStepWettingFlux(Scalar timestepFactorIn,Scalar timestepFactorOutW , Scalar& residualSatW, Scalar& residualSatNW, int globalIdxI)
{
    // compute dt restriction
    Scalar volumeCorrectionFactorIn = (1-residualSatW - residualSatNW);
    Scalar volumeCorrectionFactorOutW = (this->transProblem.variables().saturation()[globalIdxI]-residualSatW);

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (volumeCorrectionFactorOutW <= 0)
    {
        volumeCorrectionFactorOutW = 1e100;
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

    timestepFactorIn = volumeCorrectionFactorIn/timestepFactorIn;
    timestepFactorOutW = volumeCorrectionFactorOutW/timestepFactorOutW;

    Scalar timestepFactor = std::min(timestepFactorIn,timestepFactorOutW);

    return timestepFactor;
}
template<class GridView, class Scalar, class VC, class Problem> void FVSaturationWetting2P<GridView, Scalar, VC, Problem>::updateMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        int globalIdx = this->transProblem.variables().indexDiffusion(*eIt);

        Scalar sat = this->transProblem.variables().saturation()[globalIdx];

        std::vector<Scalar> mobilities = this->transProblem.materialLaw().mob(sat, globalPos, *eIt, localPos);

        // initialize mobilities
        this->transProblem.variables().mobilityWetting()[globalIdx]= mobilities[0];
        this->transProblem.variables().mobilityNonWetting()[globalIdx]= mobilities[1];
        this->transProblem.variables().capillaryPressure()[globalIdx]= this->transProblem.materialLaw().pC(sat, globalPos, *eIt, localPos);
        this->transProblem.variables().fracFlowFuncWetting()[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
        this->transProblem.variables().fracFlowFuncNonWetting()[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
    }
    return;
}

}
#endif
