// $Id: fvsaturation2p.hh 3784 2010-06-24 13:43:57Z bernd $
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
#ifndef DUMUX_FVSATURATION2P_HH
#define DUMUX_FVSATURATION2P_HH

#include "dumux/decoupled/2p/transport/fv/diffusivepart.hh"
#include "dumux/decoupled/2p/transport/fv/convectivepart.hh"
#include <dumux/decoupled/2p/transport/transportproperties.hh>
#include <dumux/decoupled/2p/2pproperties.hh>

/**
 * @file
 * @brief  Finite Volume discretization of the non-wetting phase saturation equation
 * @author Markus Wolff
 */

namespace Dumux
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
template<class TypeTag>
class FVSaturation2P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) ReferenceElements;
    typedef typename ReferenceElements::Container ReferenceElementContainer;
    typedef typename ReferenceElements::ContainerFaces ReferenceElementFaceContainer;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DiffusivePart)) DiffusivePart;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConvectivePart)) ConvectivePart;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution RepresentationType;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    DiffusivePart& diffusivePart()
    {
        return *diffusivePart_;
    }

    const DiffusivePart& diffusivePart() const
    {
        return *diffusivePart_;
    }

    ConvectivePart& convectivPart()
    {
        return *convectivePart_;
    }

    const ConvectivePart& convectivePart() const
    {
        return *convectivePart_;
    }

    //function to calculate the time step if a non-wetting phase velocity is used
    Scalar evaluateTimeStepPhaseFlux(Scalar timestepFactorIn, Scalar timestepFactorOutNW, Scalar& residualSatW,
            Scalar& residualSatNW, int globalIdxI);

    //function to calculate the time step if a total velocity is used
    Scalar evaluateTimeStepTotalFlux(Scalar timestepFactorIn, Scalar timestepFactorOut, Scalar diffFactorIn,
            Scalar diffFactorOut, Scalar& residualSatW, Scalar& residualSatNW);

public:
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
    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac, bool impes);

    //! Sets the initial solution \f$S_0\f$.
    void initialTransport();

    //! Update the values of the material laws and constitutive relations.
    /*!
     *  Constitutive relations like capillary pressure-saturation relationships, mobility-saturation relationships... are updated and stored in the variable class
     *  of type Dumux::VariableClass2P. The update has to be done when new saturation are available.
     */
    void updateMaterialLaws(RepresentationType& saturation, bool iterate);

    //! Constructs a FVSaturationNonWetting2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param velocityType a string giving the type of velocity used (could be: vn, vt)
     * \param diffPart a object of class Dune::DiffusivePart or derived from Dune::DiffusivePart (only used with vt)
     */

    FVSaturation2P(Problem& problem) :
        problem_(problem), switchNormals_(false)
    {
        if (compressibility_ && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (saturationType_ != Sw && saturationType_ != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }


        diffusivePart_ = new DiffusivePart(problem);
        convectivePart_ = new ConvectivePart(problem);
    }

    ~FVSaturation2P()
    {
        delete diffusivePart_;
        delete convectivePart_;
    }

private:
    Problem& problem_;
    DiffusivePart* diffusivePart_;
    ConvectivePart* convectivePart_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility));
    ;
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation));
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation));
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));
    bool switchNormals_;
};

template<class TypeTag>
int FVSaturation2P<TypeTag>::update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac = 1,
        bool impes = false)
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
    Dune::FieldVector<Scalar, dimWorld> gravity = problem_.gravity();

    // compute update vector
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition &localPos = ReferenceElementContainer::general(gt).position(0, 0);

        //
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        Scalar residualSatW = problem_.spatialParameters().materialLawParams(globalPos, *eIt).Swr();
        Scalar residualSatNW = problem_.spatialParameters().materialLawParams(globalPos, *eIt).Snr();

        //for benchmark only!
        //        problem_.variables().storeSrn(residualSatNW, globalIdxI);
        //for benchmark only!

        Scalar porosity = problem_.spatialParameters().porosity(globalPos, *eIt);

        Scalar viscosityWI = problem_.variables().viscosityWetting(globalIdxI);
        Scalar viscosityNWI = problem_.variables().viscosityNonwetting(globalIdxI);
        Scalar viscosityRatio = 1 - fabs(0.5 - viscosityNWI / (viscosityWI + viscosityNWI));//1 - fabs(viscosityWI-viscosityNWI)/(viscosityWI+viscosityNWI);

        Scalar densityWI = problem_.variables().densityWetting(globalIdxI);
        Scalar densityNWI = problem_.variables().densityNonwetting(globalIdxI);

        Scalar lambdaWI = problem_.variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = problem_.variables().mobilityNonwetting(globalIdxI);

        Scalar timestepFactorIn = 0;
        Scalar timestepFactorOut = 0;
        Scalar diffFactorIn = 0;
        Scalar diffFactorOut = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // local number of facet
            int isIndex = isIt->indexInInside();

            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar, dim - 1>& faceLocal =
                    ReferenceElementFaceContainer::general(faceGT).position(0, 0);

            // center of face inside volume reference element
            const LocalPosition localPosFace(0);

            Dune::FieldVector<Scalar, dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);
            if (switchNormals_)
                unitOuterNormal *= -1.0;

            Scalar faceArea = isIt->geometry().volume();

            Scalar factor = 0;
            Scalar factorSecondPhase = 0;

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                // compute factor in neighbor
                Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition& localPosNeighbor = ReferenceElementContainer::general(neighborGT).position(0, 0);

                // cell center in global coordinates
                const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                Dune::FieldVector<Scalar, dimWorld> distVec = globalPosNeighbor - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                unitDistVec /= dist;

                //get phase potentials
                Scalar potentialW = problem_.variables().potentialWetting(globalIdxI, isIndex);
                Scalar potentialNW = problem_.variables().potentialNonwetting(globalIdxI, isIndex);

                Scalar densityWJ = problem_.variables().densityWetting(globalIdxJ);
                Scalar densityNWJ = problem_.variables().densityNonwetting(globalIdxJ);

                //get velocity*normalvector*facearea/(volume*porosity)
                factor = (problem_.variables().velocity()[globalIdxI][isIndex] * unitDistVec) * (faceArea
                        * (unitOuterNormal * unitDistVec)) / (volume * porosity);
                factorSecondPhase = (problem_.variables().velocitySecondPhase()[globalIdxI][isIndex] * unitDistVec)
                        * (faceArea * (unitOuterNormal * unitDistVec)) / (volume * porosity);

                Scalar lambdaW = 0;
                Scalar lambdaNW = 0;
                Scalar viscosityW = 0;
                Scalar viscosityNW = 0;

                //upwinding of lambda dependend on the phase potential gradients
                if (potentialW >= 0.)
                {
                    lambdaW = lambdaWI;
                    if (compressibility_)
                    {
                        lambdaW /= densityWI;
                    }//divide by density because lambda is saved as lambda*density
                    viscosityW = viscosityWI;
                }
                else
                {
                    lambdaW = problem_.variables().mobilityWetting(globalIdxJ);
                    if (compressibility_)
                    {
                        lambdaW /= densityWJ;
                    }//divide by density because lambda is saved as lambda*density
                    viscosityW = problem_.variables().viscosityWetting(globalIdxJ);
                }

                if (potentialNW >= 0.)
                {
                    lambdaNW = lambdaNWI;
                    if (compressibility_)
                    {
                        lambdaNW /= densityNWI;
                    }//divide by density because lambda is saved as lambda*density
                    viscosityNW = viscosityNWI;
                }
                else
                {
                    lambdaNW = problem_.variables().mobilityNonwetting(globalIdxJ);
                    if (compressibility_)
                    {
                        lambdaNW /= densityNWJ;
                    }//divide by density because lambda is saved as lambda*density
                    viscosityNW = problem_.variables().viscosityNonwetting(globalIdxJ);
                }
                Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;

                switch (velocityType_)
                {
                case vt:
                {
                    //for time step criterion

                    if (factor >= 0)
                    {
                        timestepFactorOut += factor / (krSum * viscosityRatio);
                    }
                    if (factor < 0)
                    {
                        timestepFactorIn -= factor / (krSum * viscosityRatio);
                    }

                    //determine phase saturations from primary saturation variable
                    Scalar satWI = 0;
                    Scalar satWJ = 0;
                    switch (saturationType_)
                    {
                    case Sw:
                    {
                        satWI = problem_.variables().saturation()[globalIdxI];
                        satWJ = problem_.variables().saturation()[globalIdxJ];
                        break;
                    }
                    case Sn:
                    {
                        satWI = 1 - problem_.variables().saturation()[globalIdxI];
                        satWJ = 1 - problem_.variables().saturation()[globalIdxJ];
                        break;
                    }
                    }

                    Scalar pcI = problem_.variables().capillaryPressure(globalIdxI);
                    Scalar pcJ = problem_.variables().capillaryPressure(globalIdxJ);

                    // calculate the saturation gradient
                    Dune::FieldVector<Scalar, dimWorld> pcGradient = unitDistVec;
                    pcGradient *= (pcI - pcJ) / dist;

                    // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
                    Scalar diffPart = diffusivePart()(*eIt, isIndex, satWI, satWJ, pcGradient) * unitDistVec * faceArea
                            / (volume * porosity) * (unitOuterNormal * unitDistVec);
                    Scalar convPart = convectivePart() (*eIt, isIndex, satWI, satWJ) * unitDistVec * faceArea / (volume * porosity) * (unitOuterNormal * unitDistVec);

                    //for time step criterion
                    if (diffPart >= 0)
                    {
                        diffFactorIn += diffPart / (krSum * viscosityRatio);
                    }
                    if (diffPart < 0)
                    {
                        diffFactorOut -= diffPart / (krSum * viscosityRatio);
                    }
                    if (convPart >= 0)
                    {
                        timestepFactorOut += convPart / (krSum * viscosityRatio);
                    }
                    if (convPart < 0)
                    {
                        timestepFactorIn -= convPart / (krSum * viscosityRatio);
                    }

                    switch (saturationType_)
                    {
                    case Sw:
                    {
                        //vt*fw
                        factor *= lambdaW / (lambdaW + lambdaNW);
                        break;
                    }
                    case Sn:
                    {
                        //vt*fn
                        factor *= lambdaNW / (lambdaW + lambdaNW);
                        break;
                    }
                    }
                    factor -= diffPart;
                    factor += convPart;
                    break;
                }
                case vw:
                {
                    if (compressibility_)
                    {
                        factor /= densityWI;
                        factorSecondPhase /= densityNWI;
                    }

                    //for time step criterion
                    if (factor >= 0)
                    {
                        timestepFactorOut += factor / (krSum * viscosityRatio);
                    }
                    if (factor < 0)
                    {
                        timestepFactorIn -= factor / (krSum * viscosityRatio);
                    }
                    if (factorSecondPhase < 0)
                    {
                        timestepFactorIn -= factorSecondPhase / (krSum * viscosityRatio);
                    }

                    if (std::isnan(timestepFactorIn) || std::isinf(timestepFactorIn))
                    {
                        timestepFactorIn = 1e-100;
                    }
                    break;
                }

                    //for time step criterion if the non-wetting phase velocity is used
                case vn:
                {
                    if (compressibility_)
                    {
                        factor /= densityNWI;
                        factorSecondPhase /= densityWI;
                    }

                    //for time step criterion
                    if (factor >= 0)
                    {
                        timestepFactorOut += factor / (krSum * viscosityRatio);
                    }
                    if (factor < 0)
                    {
                        timestepFactorIn -= factor / (krSum * viscosityRatio);
                    }
                    if (factorSecondPhase < 0)
                    {
                        timestepFactorIn -= factorSecondPhase / (krSum * viscosityRatio);
                    }

                    if (std::isnan(timestepFactorIn) || std::isinf(timestepFactorIn))
                    {
                        timestepFactorIn = 1e-100;
                    }
                    break;
                }
                }
            }//end intersection with neighbor element

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                // cell center in global coordinates
                GlobalPosition globalPos = eIt->geometry().global(localPos);

                // distance vector between barycenters
                Dune::FieldVector<Scalar, dimWorld> distVec = globalPosFace - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                unitDistVec /= dist;

                //get boundary type
                BoundaryConditions::Flags bcTypeSat = problem_.bctypeSat(globalPosFace, *isIt);

                if (bcTypeSat == BoundaryConditions::dirichlet)
                {
                    Scalar satBound = problem_.dirichletSat(globalPosFace, *isIt);

                    //get velocity*normalvector*facearea/(volume*porosity)
                    factor = (problem_.variables().velocity()[globalIdxI][isIndex] * unitDistVec) * (faceArea
                            * (unitOuterNormal * unitDistVec)) / (volume * porosity);
                    factorSecondPhase = (problem_.variables().velocitySecondPhase()[globalIdxI][isIndex] * unitDistVec)
                            * (faceArea * (unitOuterNormal * unitDistVec)) / (volume * porosity);

                    Scalar pressBound = problem_.variables().pressure()[globalIdxI];
                    Scalar temperature = problem_.temperature(globalPosFace, *eIt);

                    //determine phase saturations from primary saturation variable
                    Scalar satWI = 0;
                    Scalar satWBound = 0;
                    switch (saturationType_)
                    {
                    case Sw:
                    {
                        satWI = problem_.variables().saturation()[globalIdxI];
                        satWBound = satBound;
                        break;
                    }
                    case Sn:
                    {
                        satWI = 1 - problem_.variables().saturation()[globalIdxI];
                        satWBound = 1 - satBound;
                        break;
                    }
                    }

                    Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
                            satBound);

                    //determine phase pressures from primary pressure variable
                    Scalar pressW = 0;
                    Scalar pressNW = 0;
                    switch (pressureType_)
                    {
                    case pw:
                    {
                        pressW = pressBound;
                        pressNW = pressBound + pcBound;
                        break;
                    }
                    case pn:
                    {
                        pressW = pressBound - pcBound;
                        pressNW = pressBound;
                        break;
                    }
                    }

                    FluidState fluidState;
                    fluidState.update(satBound, pressW, pressNW, temperature);

                    //get phase potentials
                    Scalar potentialW = problem_.variables().potentialWetting(globalIdxI, isIndex);
                    Scalar potentialNW = problem_.variables().potentialNonwetting(globalIdxI, isIndex);

                    Scalar lambdaW = 0;
                    Scalar lambdaNW = 0;

                    //upwinding of lambda dependend on the phase potential gradients
                    if (potentialW > 0.)
                    {
                        lambdaW = lambdaWI;
                        if (compressibility_)
                        {
                            lambdaW /= densityWI;
                        }//divide by density because lambda is saved as lambda*density
                    }
                    else
                    {
                        if (compressibility_)
                        {
                            lambdaW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
                                    satWBound) / FluidSystem::phaseViscosity(wPhaseIdx, temperature, pressW, fluidState);
                        }
                        else
                        {
                            lambdaW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt),
                                    satWBound) / viscosityWI;
                        }
                    }

                    if (potentialNW >= 0.)
                    {
                        lambdaNW = lambdaNWI;
                        if (compressibility_)
                        {
                            lambdaNW /= densityNWI;
                        }//divide by density because lambda is saved as lambda*density
                    }
                    else
                    {
                        if (compressibility_)
                        {
                            lambdaNW = MaterialLaw::krn(
                                    problem_.spatialParameters().materialLawParams(globalPos, *eIt), satWBound)
                                    / FluidSystem::phaseViscosity(nPhaseIdx, temperature, pressNW, fluidState);
                        }
                        else
                        {
                            lambdaNW = MaterialLaw::krn(
                                    problem_.spatialParameters().materialLawParams(globalPos, *eIt), satWBound)
                                    / viscosityNWI;
                        }
                    }
                    //                    std::cout<<lambdaW<<" "<<lambdaNW<<std::endl;

                    Scalar krSum = lambdaW * viscosityWI + lambdaNW * viscosityNWI;

                    switch (velocityType_)
                    {
                    case vt:
                    {
                        //for time step criterion

                        if (factor >= 0)
                        {
                            timestepFactorOut += factor / (krSum * viscosityRatio);
                        }
                        if (factor < 0)
                        {
                            timestepFactorIn -= factor / (krSum * viscosityRatio);
                        }

                        Scalar pcI = problem_.variables().capillaryPressure(globalIdxI);

                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar, dimWorld> pcGradient = unitDistVec;
                        pcGradient *= (pcI - pcBound) / dist;

                        // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
                        Scalar diffPart = diffusivePart()(*eIt, isIndex, satWI, satWBound, pcGradient) * unitDistVec
                                * faceArea / (volume * porosity) * (unitOuterNormal * unitDistVec);
                        Scalar convPart = convectivePart()(*eIt, isIndex, satWI, satWBound) * unitDistVec * faceArea / (volume * porosity) * (unitOuterNormal * unitDistVec);

                        //for time step criterion
                        if (diffPart >= 0)
                        {
                            diffFactorIn += diffPart / (krSum * viscosityRatio);
                        }
                        if (diffPart < 0)
                        {
                            diffFactorOut -= diffPart / (krSum * viscosityRatio);
                        }
                        if (convPart >= 0)
                        {
                            timestepFactorOut += convPart / (krSum * viscosityRatio);
                        }
                        if (convPart < 0)
                        {
                            timestepFactorIn -= convPart / (krSum * viscosityRatio);
                        }

                        switch (saturationType_)
                        {
                        case Sw:
                        {
                            //vt*fw
                            factor *= lambdaW / (lambdaW + lambdaNW);
                            break;
                        }
                        case Sn:
                        {
                            //vt*fn
                            factor *= lambdaNW / (lambdaW + lambdaNW);
                            break;
                        }
                        }
                        //vt*fw
                        factor -= diffPart;
                        factor += convPart;
                        break;
                    }

                    case vw:
                    {
                        if (compressibility_)
                        {
                            factor /= densityWI;
                            factorSecondPhase /= densityNWI;
                        }

                        //for time step criterion
                        if (factor >= 0)
                        {
                            timestepFactorOut += factor / (krSum * viscosityRatio);
                        }
                        if (factor < 0)
                        {
                            timestepFactorIn -= factor / (krSum * viscosityRatio);
                        }
                        if (factorSecondPhase < 0)
                        {
                            timestepFactorIn -= factorSecondPhase / (krSum * viscosityRatio);
                        }

                        if (std::isnan(timestepFactorIn) || std::isinf(timestepFactorIn))
                        {
                            timestepFactorIn = 1e-100;
                        }
                        break;
                    }

                        //for time step criterion if the non-wetting phase velocity is used
                    case vn:
                    {
                        if (compressibility_)
                        {
                            factor /= densityNWI;
                            factorSecondPhase /= densityWI;
                        }

                        //for time step criterion
                        if (factor >= 0)
                        {
                            timestepFactorOut += factor / (krSum * viscosityRatio);
                        }
                        if (factor < 0)
                        {
                            timestepFactorIn -= factor / (krSum * viscosityRatio);
                        }
                        if (factorSecondPhase < 0)
                        {
                            timestepFactorIn -= factorSecondPhase / (krSum * viscosityRatio);
                        }

                        if (std::isnan(timestepFactorIn) || std::isinf(timestepFactorIn))
                        {
                            timestepFactorIn = 1e-100;
                        }
                        break;
                    }
                    }
                }//end dirichlet boundary

                if (bcTypeSat == BoundaryConditions::neumann)
                {
                    //get mobilities
                    Scalar lambdaW, lambdaNW;

                    lambdaW = lambdaWI;
                    if (compressibility_)
                    {
                        lambdaW /= densityWI;
                    }

                    lambdaNW = lambdaNWI;
                    if (compressibility_)
                    {
                        lambdaNW /= densityNWI;
                    }

                    Scalar krSum = lambdaW * viscosityWI + lambdaNW * viscosityNWI;

                    //get velocity*normalvector*facearea/(volume*porosity)
                    factor = (problem_.variables().velocity()[globalIdxI][isIndex] * unitOuterNormal) * faceArea
                            / (volume * porosity);
                    switch (velocityType_)
                    {
                    case vt:
                    {
                        switch (saturationType_)
                        {
                        case Sw:
                        {
                            //vt*fw
                            factor *= lambdaW / (lambdaW + lambdaNW);
                            break;
                        }
                        case Sn:
                        {
                            //vt*fn
                            factor *= lambdaNW / (lambdaW + lambdaNW);
                            break;
                        }
                        }
                        break;
                    }
                    }
                    Scalar boundaryFactor = problem_.neumannSat(globalPosFace, *isIt, factor);

                    if (factor != boundaryFactor)
                    {
                        switch (saturationType_)
                        {
                        case Sw:
                        {
                            factor = boundaryFactor / densityWI * faceArea / (volume * porosity);
                            break;
                        }
                        case Sn:
                        {
                            factor = boundaryFactor / densityNWI * faceArea / (volume * porosity);
                            break;
                        }
                        }
                    }

                    //for time step criterion

                    if (factor >= 0)
                    {
                        timestepFactorOut += factor / (krSum * viscosityRatio);
                    }
                    if (factor < 0)
                    {
                        timestepFactorIn -= factor / (krSum * viscosityRatio);
                    }
                }//end neumann boundary
            }//end boundary
            // add to update vector
            updateVec[globalIdxI] -= factor;//-:v>0, if flow leaves the cell
        }// end all intersections
        Scalar source = 0;
        switch (velocityType_)
        {
        case vw:
        {
            source = problem_.source(globalPos, *eIt)[wPhaseIdx] / densityWI;
            break;
        }
        case vn:
        {
            source = problem_.source(globalPos, *eIt)[nPhaseIdx] / densityNWI;
            break;
        }
        case vt:
        {
            source = problem_.source(globalPos, *eIt)[wPhaseIdx] / densityWI + problem_.source(globalPos,
                    *eIt)[nPhaseIdx] / densityNWI;
            break;
        }
        }
        if (source)
        {
            //get mobilities
            Scalar lambdaW = 0;
            Scalar lambdaNW = 0;

            lambdaW = lambdaWI;
            if (compressibility_)
            {
                lambdaW /= densityWI;
            }
            lambdaNW = lambdaNWI;
            if (compressibility_)
            {
                lambdaNW /= densityNWI;
            }

            Scalar krSum = lambdaW * viscosityWI + lambdaNW * viscosityNWI;
            switch (saturationType_)
            {
            case Sw:
            {
                updateVec[globalIdxI] += source / porosity;
                break;
            }
            case Sn:
            {
                updateVec[globalIdxI] += source / porosity;
                break;
            }
            }
            if (source >= 0)
            {
                timestepFactorIn += source / (porosity * viscosityRatio * krSum);
            }
            else
            {
                timestepFactorOut -= source / (porosity * viscosityRatio * krSum);
            }
        }

        //calculate time step
        if (velocityType_ == vw || velocityType_ == vn)
        {
            dt = std::min(dt, evaluateTimeStepPhaseFlux(timestepFactorIn, timestepFactorOut, residualSatW,
                    residualSatNW, globalIdxI));
        }
        if (velocityType_ == vt)
        {
            dt = std::min(dt, evaluateTimeStepTotalFlux(timestepFactorIn, timestepFactorOut, diffFactorIn,
                    diffFactorOut, residualSatW, residualSatNW));
        }

        problem_.variables().volumecorrection(globalIdxI) = updateVec[globalIdxI];
    } // end grid traversal

    return 0;
}

template<class TypeTag>
void FVSaturation2P<TypeTag>::initialTransport()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition &localPos = ReferenceElementContainer::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // initialize cell concentration
        problem_.variables().saturation()[problem_.variables().index(*eIt)] = problem_.initSat(globalPos, *eIt);
    }
    return;
}

template<class TypeTag>
typename FVSaturation2P<TypeTag>::Scalar FVSaturation2P<TypeTag>::evaluateTimeStepTotalFlux(Scalar timestepFactorIn,
        Scalar timestepFactorOut, Scalar diffFactorIn, Scalar diffFactorOut, Scalar& residualSatW,
        Scalar& residualSatNW)
{
    // compute volume correction
    Scalar volumeCorrectionFactor = (1 - residualSatW - residualSatNW);

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (timestepFactorIn <= 0)
    {
        timestepFactorIn = 1e-100;
    }
    if (timestepFactorOut <= 0)
    {
        timestepFactorOut = 1e-100;
    }

    Scalar sumFactor = std::min(volumeCorrectionFactor / timestepFactorIn, volumeCorrectionFactor / timestepFactorOut);

    //make sure that diffFactor > 0
    if (diffFactorIn <= 0)
    {
        diffFactorIn = 1e-100;
    }
    if (diffFactorOut <= 0)
    {
        diffFactorOut = 1e-100;
    }

    Scalar minDiff = std::min(volumeCorrectionFactor / diffFactorIn, volumeCorrectionFactor / diffFactorOut);

    //determine time step
    sumFactor = std::min(sumFactor, 0.1 * minDiff);

    return sumFactor;
}
template<class TypeTag>
typename FVSaturation2P<TypeTag>::Scalar FVSaturation2P<TypeTag>::evaluateTimeStepPhaseFlux(
        Scalar timestepFactorIn, Scalar timestepFactorOut, Scalar& residualSatW, Scalar& residualSatNW, int globalIdxI)
{
    // compute dt restriction
    Scalar volumeCorrectionFactorIn = 1 - residualSatW - residualSatNW;
    Scalar volumeCorrectionFactorOut = 0;
    if (saturationType_ == Sw)
    {
        Scalar satI = problem_.variables().saturation()[globalIdxI];
        volumeCorrectionFactorOut = std::max((satI - residualSatW), 1e-2);
    }
    if (saturationType_ == Sn)
    {
        Scalar satI = problem_.variables().saturation()[globalIdxI];
        volumeCorrectionFactorOut = std::max((satI - residualSatNW), 1e-2);
    }

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
    if (timestepFactorOut <= 0)
    {
        timestepFactorOut = 1e-100;
    }

    //correct volume
    timestepFactorIn = volumeCorrectionFactorIn / timestepFactorIn;
    timestepFactorOut = volumeCorrectionFactorOut / timestepFactorOut;

    //determine timestep
    Scalar timestepFactor = std::min(timestepFactorIn, timestepFactorOut);

    return timestepFactor;
}
template<class TypeTag>
void FVSaturation2P<TypeTag>::updateMaterialLaws(RepresentationType& saturation = *(new RepresentationType(0)),
        bool iterate = false)
{
    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition &localPos = ReferenceElementContainer::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        int globalIdx = problem_.variables().index(*eIt);

        Scalar sat = 0;
        if (!iterate)
        {
            sat = problem_.variables().saturation()[globalIdx];
        }
        else
        {
            sat = saturation[globalIdx];
        }

        //determine phase saturations from primary saturation variable
        Scalar satW = 0;
        if (saturationType_ == Sw)
        {
            satW = sat;
        }
        if (saturationType_ == Sn)
        {
            satW = 1 - sat;
        }

        Scalar temperature = problem_.temperature(globalPos, *eIt);
        Scalar referencePressure =  problem_.referencePressure(globalPos, *eIt);

        fluidState.update(satW, referencePressure, referencePressure, temperature);

        problem_.variables().densityWetting(globalIdx) = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
        problem_.variables().densityNonwetting(globalIdx) = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
        problem_.variables().viscosityWetting(globalIdx) = FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
        problem_.variables().viscosityNonwetting(globalIdx) = FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx) = MaterialLaw::krw(
                problem_.spatialParameters().materialLawParams(globalPos, *eIt), satW)
                / problem_.variables().viscosityWetting(globalIdx);
        problem_.variables().mobilityNonwetting(globalIdx) = MaterialLaw::krn(
                problem_.spatialParameters().materialLawParams(globalPos, *eIt), satW)
                / problem_.variables().viscosityNonwetting(globalIdx);
        problem_.variables().capillaryPressure(globalIdx) = MaterialLaw::pC(
                problem_.spatialParameters().materialLawParams(globalPos, *eIt), satW);

        problem_.variables().fracFlowFuncWetting(globalIdx) = problem_.variables().densityWetting(globalIdx)
                * problem_.variables().mobilityWetting(globalIdx) / (problem_.variables().densityWetting(globalIdx)
                * problem_.variables().mobilityWetting(globalIdx) + problem_.variables().densityNonwetting(globalIdx)
                * problem_.variables().mobilityNonwetting(globalIdx));
        problem_.variables().fracFlowFuncNonwetting(globalIdx) = problem_.variables().densityNonwetting(globalIdx)
                * problem_.variables().mobilityNonwetting(globalIdx) / (problem_.variables().densityWetting(globalIdx)
                * problem_.variables().mobilityWetting(globalIdx) + problem_.variables().densityNonwetting(globalIdx)
                * problem_.variables().mobilityNonwetting(globalIdx));

        problem_.spatialParameters().update(satW, *eIt);
    }
    return;
}

}
#endif
