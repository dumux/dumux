// $Id:$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_FVSATURATIONNONWETTING2P_HH
#define DUNE_FVSATURATIONNONWETTING2P_HH

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
 * @brief  Finite Volume discretization of the non-wetting phase saturation equation
 * @author Markus Wolff
 */

namespace Dune
{
//! \ingroup transport
//! The finite volume model for the solution of the non-wetting phase saturation equation
/*! Provides a Finite Volume implementation for the evaluation
 *  of equations of the form
 *  \f[
 *    \frac{\partial S_n}{\partial t} + \text{div}\, \boldsymbol{v_n} = 0,
 *  \f]
 *  where \f$\boldsymbol{v}_n = \lambda_n \boldsymbol{K} \left(\text{grad}\, p_n + \rho_n g  \text{grad}\, z\right)\f$,
 *  where \f$p_n\f$ denotes the wetting phase pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda_n\f$ the non-wetting phase mobility,
 *  \f$\rho_n\f$ the non-wetting phase density and \f$g\f$ the gravity constant and \f$S_n\f$ the non-wetting phase saturation,
 *
 *  or where \f$\boldsymbol{v}_n = f_n \boldsymbol{v_{total}} - f_n \lambda_w \boldsymbol{K} \text{grad}\, p_c \f$,
 *  \f$f_n\f$ is the non-wetting phase fractional flow function, \f$\lambda_w\f$ is the wetting phase mobility, \f$\boldsymbol{K}\f$ the absolute permeability,
 *  \f$p_c\f$ the capillary pressure and \f$S_n\f$ the wetting phase saturation.
 *
 *

 Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem

 */
template<class GridView, class Scalar, class VC,
        class Problem = TransportProblem<GridView, Scalar, VC> >
class FVSaturationNonWetting2P: public Transport<GridView, Scalar, VC, Problem>
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
        vn = 0, vt = 1
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

    //function to calculate the time step if a non-wetting phase velocity is used
    Scalar evaluateTimeStepNonWettingFlux(Scalar timestepFactorIn, Scalar timestepFactorOutNW, Scalar& residualSatW, Scalar& residualSatNW, int globalIdxI);

    //function to calculate the time step if a total velocity is used
    Scalar evaluateTimeStepTotalFlux(Scalar timestepFactorIn,Scalar timestepFactorOut, Scalar diffFactorIn, Scalar diffFactorOut, Scalar& residualSatW, Scalar& residualSatNW);

public:
    typedef typename VC::ScalarVectorType RepresentationType;//!< Data type for a Vector of Scalars

    //! Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[in] dt         time step size
     *  \param[in] updateVec  vector for the update values
     *  \param[in] CLFFac     security factor for the time step criterion (0 < CLFFac <= 1)
     *  \param[in] impes      variable is true if an impes algorithm is used and false if the transport part is solved independently
     *
     *  This method calculates the update vector \f$ u \f$ of the discretized equation
     *  \f[
     *   S_{n_{new}} = S_{n_{old}} - u,
     *  \f]
     *  where \f$ u = \sum_{element faces} \boldsymbol{v}_n * \boldsymbol{n} * A_{element face}\f$, \f$\boldsymbol{n}\f$ is the face normal and \f$A_{element face}\f$ is the face area.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing a CFL condition.
     */
    int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac, bool impes);

    //! Sets the initial solution \f$S_0\f$.
    void initialTransport();

    //! Update the values of the material laws and constitutive relations.
    /*!
     *  Constitutive relations like capillary pressure-saturation relationships, mobility-saturation relationships... are updated and stored in the variable class
     *  of type Dune::VariableClass2P. The update has to be done when new saturation are available.
     */
    void updateMaterialLaws(RepresentationType& saturation, bool iterate);

    //! Write data files
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    virtual void vtkout(const char* name, int k) const
    {
        this->transProblem.variables().vtkout(name, k);
        return;
    }

    //! Constructs a FVSaturationNonWetting2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param velocityType a string giving the type of velocity used (could be: vn, vt)
     * \param diffPart a object of class Dune::DiffusivePart or derived from Dune::DiffusivePart (only used with vt)
     */
    FVSaturationNonWetting2P(GridView& gridView, Problem& problem, std::string velocityType, DiffusivePart<GridView,Scalar>& diffPart = *(new DiffusivePart<GridView, Scalar>))
    :Transport<GridView, Scalar, VC, Problem>(gridView, problem), diffusivePart_(diffPart),
    velocityType_((velocityType == "vn") ? 0 : ((velocityType == "vt") ? 1 : 999))
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
int FVSaturationNonWetting2P<GridView, Scalar, VC, Problem>::update(const Scalar t, Scalar& dt,
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

                //get phase potentials
                Scalar potentialW = this->transProblem.variables().potentialWetting()[globalIdxI][indexInInside];
                Scalar potentialNW = this->transProblem.variables().potentialNonWetting()[globalIdxI][indexInInside];

                //get velocity*normalvector*facearea/(volume*porosity)
                factor = (this->transProblem.variables().velocity()[globalIdxI][indexInInside] * unitOuterNormal) * faceArea / (volume*porosity);

                if (velocityType_ == vt)
                {
                    Scalar lambdaW, lambdaNW;

                    //upwinding of lambda dependend on the phase potential gradients
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

                    // calculate the saturation gradient
                    Dune::FieldVector<Scalar,dimWorld> pcGradient = distVec;
                    pcGradient *= (pcJ - pcI)/(dist*dist);

                    // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
                    Scalar diffPart = diffusivePart_(*eIt, indexInInside, 1-satI, 1-satJ, pcGradient)*unitOuterNormal * faceArea / (volume*porosity);

                    //for time step criterion
                    if (diffPart >= 0)
                    {
                        diffFactorIn += diffPart/(krSum*viscosityRatio);
                    }
                    if (diffPart < 0)
                    {
                        diffFactorOut -= diffPart/(krSum*viscosityRatio);
                    }

                    //vt*fw
                    factor *= lambdaNW/(lambdaW + lambdaNW);
                    factor -= diffPart;
                }

                //for time step criterion if the non-wetting phase velocity is used
                if (velocityType_ == vn)
                {
                    if (potentialNW >= 0)
                    {
                        timestepFactorOut+= factor;

                        if (potentialW < 0)
                        {
                            Scalar lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxJ];
                            Scalar lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                            Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                            //get guess of wetting flux by weighting of the non-wetting flux with the potential gradients
                            timestepFactorIn -= factor/(lambdaNW*potentialNW*viscosityRatio*krSum)*lambdaW*potentialW;
                        }

                    }
                    if (potentialNW < 0)
                    {
                        Scalar lambdaW = 0;
                        Scalar lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxJ];
                        Scalar krSum = 1;
                        if (potentialW >= 0)
                        {
                            lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                            krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                        }
                        if (potentialW < 0)
                        {
                            lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxJ];
                            krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                            //get guess of wetting flux by weighting of the non-wetting flux with the potential gradients
                            timestepFactorIn -= factor/(lambdaNW*potentialNW*viscosityRatio*krSum)*lambdaW*potentialW;
                        }

                        timestepFactorIn -= factor/(viscosityRatio*krSum);

                    }
                    if (std::isnan(timestepFactorIn) || std::isinf(timestepFactorIn))
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

                    //get phase potentials
                    Scalar potentialW = this->transProblem.variables().potentialWetting()[globalIdxI][indexInInside];
                    Scalar potentialNW = this->transProblem.variables().potentialNonWetting()[globalIdxI][indexInInside];

                    //get velocity*normalvector*facearea/(volume*porosity)
                    factor = (this->transProblem.variables().velocity()[globalIdxI][indexInInside] * unitOuterNormal) * faceArea / (volume*porosity);

                    if (velocityType_ == vt)
                    {
                        Scalar lambdaW, lambdaNW;

                        //upwinding of lambda dependend on the phase potential gradients
                        if (potentialW >= 0.)
                        {
                            lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                        }
                        else
                        {
                            //sat = S_n and lambda = lambda(S_W)!
                            lambdaW = this->transProblem.materialLaw().mobW(1-satBound,globalPosFace, *eIt, localPosFace);
                        }

                        if (potentialNW >= 0.)
                        {
                            lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                        }
                        else
                        {
                            lambdaNW = this->transProblem.materialLaw().mobN(satBound,globalPosFace, *eIt, localPosFace);
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
                        //sat = S_n and pc = pc(S_W)!
                        Scalar pcBound = this->transProblem.materialLaw().pC(1-satBound,globalPosFace, *eIt, localPosFace);
                        //                        std::cout<<"satBound = "<<satBound<<"pcBound = "<<pcBound<<"pcI = "<<pcI<<std::endl;

                        // compute distance between cell centers
                        Scalar dist = distVec.two_norm();

                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar,dimWorld> pcGradient = distVec;
                        pcGradient *= (pcBound - pcI)/(dist*dist);

                        // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
                        Scalar diffPart = diffusivePart_(*eIt, indexInInside, 1-satI, 1-satBound, pcGradient)* unitOuterNormal * faceArea / (volume*porosity);

                        //for time step criterion
                        if (diffPart >= 0)
                        {
                            diffFactorIn += diffPart/(krSum*viscosityRatio);
                        }
                        if (diffPart < 0)
                        {
                            diffFactorOut -= diffPart/(krSum*viscosityRatio);
                        }

                        //vt*fw
                        factor *= lambdaNW/(lambdaW + lambdaNW);
                        factor -= diffPart;
                    }

                    //for time step criterion if the non-wetting phase velocity is used
                    if (velocityType_ == vn)
                    {
                        if (potentialNW >= 0)
                        {
                            timestepFactorOut+= factor;

                            if (potentialW < 0)
                            {
                                Scalar lambdaW = this->transProblem.materialLaw().mobW(1-satBound,globalPosFace, *eIt, localPosFace);
                                Scalar lambdaNW = this->transProblem.variables().mobilityNonWetting()[globalIdxI];
                                Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                                //get guess of wetting flux by weighting of the non-wetting flux with the potential gradients
                                timestepFactorIn -= factor/(lambdaNW*potentialNW*viscosityRatio*krSum)*lambdaW*potentialW;
                            }

                        }
                        if (potentialNW < 0)
                        {
                            Scalar lambdaW = 0;
                            Scalar lambdaNW = this->transProblem.materialLaw().mobN(satBound,globalPosFace, *eIt, localPosFace);
                            Scalar krSum = 1;
                            if (potentialW >= 0)
                            {
                                lambdaW = this->transProblem.variables().mobilityWetting()[globalIdxI];
                                krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
                            }
                            if (potentialW < 0)
                            {
                                lambdaW = this->transProblem.materialLaw().mobW(1-satBound, globalPosFace, *eIt, localPosFace);
                                krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                                //get guess of wetting flux by weighting of the non-wetting flux with the potential gradients
                                timestepFactorIn -= factor/(lambdaNW*potentialNW*viscosityRatio*krSum)*lambdaW*potentialW;
                            }

                            timestepFactorIn -= factor/(viscosityRatio*krSum);
                        }
                        if (std::isnan(timestepFactorIn) || std::isinf(timestepFactorIn))
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

        //calculate time step
        if (velocityType_ == vn)
        {
            dt = std::min(dt, evaluateTimeStepNonWettingFlux(timestepFactorIn, timestepFactorOut, residualSatW, residualSatNW, globalIdxI));
        }
        if (velocityType_ == vt)
        {
            dt = std::min(dt, evaluateTimeStepTotalFlux(timestepFactorIn, timestepFactorOut, diffFactorIn, diffFactorOut, residualSatW, residualSatNW));
        }
    } // end grid traversal

    return 0;
}

template<class GridView, class Scalar, class VC, class Problem>
void FVSaturationNonWetting2P<GridView, Scalar, VC, Problem>::initialTransport()
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
Scalar FVSaturationNonWetting2P<GridView, Scalar, VC, Problem>::evaluateTimeStepTotalFlux(Scalar timestepFactorIn,Scalar timestepFactorOut, Scalar diffFactorIn, Scalar diffFactorOut, Scalar& residualSatW, Scalar& residualSatNW)
{
    // compute volume correction
    Scalar volumeCorrectionFactor = (1 - residualSatW -residualSatNW);

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (timestepFactorIn <= 0)
    {
        timestepFactorIn = 1e-100;
    }
    if (timestepFactorOut <= 0)
    {
        timestepFactorOut = 1e-100;
    }

    Scalar sumFactor = std::min(volumeCorrectionFactor/timestepFactorIn, volumeCorrectionFactor/timestepFactorOut);

    //make sure that diffFactor > 0
    if (diffFactorIn <= 0)
    {
        diffFactorIn = 1e-100;
    }
    if (diffFactorOut <= 0)
    {
        diffFactorOut = 1e-100;
    }

    Scalar minDiff = std::min(volumeCorrectionFactor/diffFactorIn,volumeCorrectionFactor/diffFactorOut);

    //determine time step
    sumFactor = std::min(sumFactor, 0.1*minDiff);

    return sumFactor;
}
template<class GridView, class Scalar, class VC, class Problem>
Scalar FVSaturationNonWetting2P<GridView, Scalar, VC, Problem>::evaluateTimeStepNonWettingFlux(Scalar timestepFactorIn,Scalar timestepFactorOut, Scalar& residualSatW, Scalar& residualSatNW, int globalIdxI)
{
    // compute dt restriction
    Scalar volumeCorrectionFactorIn = (1-residualSatW - residualSatNW);
    Scalar volumeCorrectionFactorOut = (this->transProblem.variables().saturation()[globalIdxI] - residualSatNW);

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (volumeCorrectionFactorOut <= 0)
    {
        volumeCorrectionFactorOut = 1e100;
    }

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (timestepFactorIn <= 0)
    {
        timestepFactorIn = 1e-100;
    }
    if (timestepFactorOut<= 0)
    {
        timestepFactorOut= 1e-100;
    }

    //correct volume
    timestepFactorIn = volumeCorrectionFactorIn / timestepFactorIn;
    timestepFactorOut= volumeCorrectionFactorOut / timestepFactorOut;

    //determine timestep
    Scalar timestepFactor = std::min(timestepFactorIn,timestepFactorOut);

    return timestepFactor;
}
template<class GridView, class Scalar, class VC, class Problem> void FVSaturationNonWetting2P<GridView, Scalar, VC, Problem>::updateMaterialLaws(RepresentationType& saturation = *(new RepresentationType(0)), bool iterate=false)
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

        Scalar sat = 0;
        if (!iterate)
        {
            sat = this->transProblem.variables().saturation()[globalIdx];
        }
        else
        {
            sat = saturation[globalIdx];
        }

        std::vector<Scalar> mobilities = this->transProblem.materialLaw().mob(1-sat, globalPos, *eIt, localPos);

        // initialize mobilities
        this->transProblem.variables().mobilityWetting()[globalIdx]= mobilities[0];
        this->transProblem.variables().mobilityNonWetting()[globalIdx]= mobilities[1];
        this->transProblem.variables().capillaryPressure()[globalIdx]= this->transProblem.materialLaw().pC(1-sat, globalPos, *eIt, localPos);
        this->transProblem.variables().fracFlowFuncWetting()[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
        this->transProblem.variables().fracFlowFuncNonWetting()[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
    }
    return;
}

}
#endif
