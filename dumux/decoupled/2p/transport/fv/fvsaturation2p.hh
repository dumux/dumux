// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVSATURATION2P_HH
#define DUMUX_FVSATURATION2P_HH

#include "dumux/decoupled/2p/transport/fv/diffusivepart.hh"
#include "dumux/decoupled/2p/transport/fv/convectivepart.hh"
#include <dumux/decoupled/2p/transport/transportproperties.hh>
#include <dumux/decoupled/2p/2pproperties.hh>
#include "evalcflflux_default.hh"

/**
 * @file
 * @brief  Finite Volume discretization of a saturation transport equation
 * @author Markus Wolff
 */

namespace Dumux
{
//! \ingroup Saturation2p
//! \brief The finite volume discretization of a saturation transport equation
/*! Provides a Finite Volume implementation for the evaluation
 *  of equations of the form
 *
 *  \f[\frac{\partial S_\alpha}{\partial t} + \text{div}\, \boldsymbol{v_\alpha} = q_\alpha,\f]
 *
 *  where \f$S_\alpha\f$ is the saturation of phase alpha (wetting (w), non-wetting (n)) and \f$\boldsymbol{v}_\alpha\f$ is the phase velocity calculated by the multi-phase Darcy equation,
 *  and of the form
 *
 * \f[\frac{\partial S_w}{\partial t} + f_w \text{div}\, \boldsymbol{v}_{t} + f_w \lambda_n \boldsymbol{K}\left(\text{grad}\, p_c + (\rho_n-\rho_w) g \text{grad} z \right)= q_\alpha,\f]
 *
 * \f[\frac{\partial S_n}{\partial t} + f_n \text{div}\, \boldsymbol{v}_{t} - f_n \lambda_w \boldsymbol{K}\left(\text{grad}\, p_c + (\rho_n-\rho_w) g \text{grad} z \right)= q_\alpha,\f]
 *
 *  where \f$f_\alpha\f$ is the fractional flow function, \f$\lambda_\alpha\f$ is the mobility, \f$\boldsymbol{K}\f$ the absolute permeability,
 *  \f$p_c\f$ the capillary pressure, \f$\rho\f$ the fluid density, \f$g\f$ the gravity constant, and \f$q\f$ the source term.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVSaturation2P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef Dune::GenericReferenceElements<Scalar, dim>
            ReferenceElementContainer;
    typedef Dune::GenericReferenceElements<Scalar, dim - 1>
            ReferenceElementFaceContainer;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DiffusivePart)) DiffusivePart;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConvectivePart))
            ConvectivePart;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(EvalCflFluxFunction))
            EvalCflFluxFunction;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters))
            SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

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

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    DiffusivePart& diffusivePart()
    {
        return *diffusivePart_;
    }

    const DiffusivePart& diffusivePart() const
    {
        return *diffusivePart_;
    }

    ConvectivePart& convectivePart()
    {
        return *convectivePart_;
    }

    const ConvectivePart& convectivePart() const
    {
        return *convectivePart_;
    }

    EvalCflFluxFunction& evalCflFluxFunction()
    {
        return *evalCflFluxFunction_;
    }

    const EvalCflFluxFunction& evalCflFluxFunction() const
    {
        return *evalCflFluxFunction_;
    }

public:
    //! Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[in] dt         time step size
     *  \param[in] updateVec  vector for the update values
     *  \param[in] impes      variable is true if an impes algorithm is used and false if the transport part is solved independently
     *
     *  This method calculates the update vector \f$\boldsymbol{u}\f$ of the discretized equation
     *
     *  \f[S_{new} = S_{old} + u,\f]
     *
     *  where \f$ u = \sum_{element faces} \boldsymbol{F}\f$ and \f$\boldsymbol{F}\f$ is the flux over an element face.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing a CFL condition.
     */
    virtual int update(const Scalar t, Scalar& dt,
            RepresentationType& updateVec, bool impes);

    //! Sets the initial solution \f$S_0\f$.
    void initialize();

    //! Update the values of the material laws and constitutive relations.
    /*!
     *  Constitutive relations like capillary pressure-saturation relationships, mobility-saturation relationships... are updated and stored in the variable class
     *  of type Dumux::VariableClass2P. The update has to be done when new saturation are available.
     */
    void updateMaterialLaws(RepresentationType& saturation, bool iterate);

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        typename Variables::ScalarSolutionType *saturation =
                writer.template createField<Scalar, 1> (
                        problem_.gridView().size(0));

        *saturation = problem_.variables().saturation();

        if (saturationType_ == Sw)
        {
            writer.addCellData(saturation, "wetting saturation");
        }
        if (saturationType_ == Sn)
        {
            writer.addCellData(saturation, "nonwetting saturation");
        }

        return;
    }

    // serialization methods
    //! Function needed for restart option.
    template<class Restarter>
    void serialize(Restarter &res)
    {
        problem_.variables().serialize<Restarter> (res);
    }

    //! Function needed for restart option.
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        problem_.variables().deserialize<Restarter> (res);
    }

    //! Constructs a FVSaturation2P object
    /**

     * \param problem a problem class object
     */

    FVSaturation2P(Problem& problem) :
        problem_(problem), switchNormals_(false)
    {
        if (compressibility_ && velocityType_ == vt)
        {
            DUNE_THROW(
                    Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (saturationType_ != Sw && saturationType_ != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (pressureType_ != pw && pressureType_ != pn && pressureType_
                != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        diffusivePart_ = new DiffusivePart(problem);
        convectivePart_ = new ConvectivePart(problem);
        evalCflFluxFunction_ = new EvalCflFluxFunction(problem);
    }

    ~FVSaturation2P()
    {
        delete diffusivePart_;
        delete convectivePart_;
        delete evalCflFluxFunction_;
    }

private:
    Problem& problem_;
    DiffusivePart* diffusivePart_;
    ConvectivePart* convectivePart_;
    EvalCflFluxFunction* evalCflFluxFunction_;

    static const bool compressibility_ =
            GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility));
    ;
    static const int saturationType_ =
            GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation));
    static const int velocityType_ =
            GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation));
    static const int pressureType_ =
            GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));
    bool switchNormals_;
};

template<class TypeTag>
int FVSaturation2P<TypeTag>::update(const Scalar t, Scalar& dt,
        RepresentationType& updateVec, bool impes = false)
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
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt
            != eItEnd; ++eIt)
    {
        //
        GlobalPosition globalPos = eIt->geometry().center();

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        //for benchmark only!
        //        problem_.variables().storeSrn(residualSatNW, globalIdxI);
        //for benchmark only!

        Scalar porosity =
                problem_.spatialParameters().porosity(globalPos, *eIt);

        Scalar viscosityWI = problem_.variables().viscosityWetting(globalIdxI);
        Scalar viscosityNWI = problem_.variables().viscosityNonwetting(
                globalIdxI);

        Scalar densityWI = problem_.variables().densityWetting(globalIdxI);
        Scalar densityNWI = problem_.variables().densityNonwetting(globalIdxI);

        Scalar lambdaWI = problem_.variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = problem_.variables().mobilityNonwetting(globalIdxI);

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt
                != isItEnd; ++isIt)
        {
            // local number of facet
            int isIndex = isIt->indexInInside();

            Dune::FieldVector<Scalar, dimWorld> unitOuterNormal =
                    isIt->centerUnitOuterNormal();
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

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor =
                        neighborPointer->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector<Scalar, dimWorld> distVec = globalPosNeighbor
                        - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                //get phase potentials
                Scalar potentialW = problem_.variables().potentialWetting(
                        globalIdxI, isIndex);
                Scalar potentialNW = problem_.variables().potentialNonwetting(
                        globalIdxI, isIndex);

                Scalar densityWJ = problem_.variables().densityWetting(
                        globalIdxJ);
                Scalar densityNWJ = problem_.variables().densityNonwetting(
                        globalIdxJ);

                //get velocity*normalvector*facearea/(volume*porosity)
                factor = (problem_.variables().velocity()[globalIdxI][isIndex]
                        * unitOuterNormal) * (faceArea);
                factorSecondPhase
                        = (problem_.variables().velocitySecondPhase()[globalIdxI][isIndex]
                                * unitOuterNormal) * (faceArea);

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
                    viscosityW = problem_.variables().viscosityWetting(
                            globalIdxJ);
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
                    lambdaNW = problem_.variables().mobilityNonwetting(
                            globalIdxJ);
                    if (compressibility_)
                    {
                        lambdaNW /= densityNWJ;
                    }//divide by density because lambda is saved as lambda*density
                    viscosityNW = problem_.variables().viscosityNonwetting(
                            globalIdxJ);
                }

                switch (velocityType_)
                {
                case vt:
                {
                    //add cflFlux for time-stepping
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, factor, *isIt);

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
                        satWI = 1
                                - problem_.variables().saturation()[globalIdxI];
                        satWJ = 1
                                - problem_.variables().saturation()[globalIdxJ];
                        break;
                    }
                    }

                    Scalar pcI = problem_.variables().capillaryPressure(
                            globalIdxI);
                    Scalar pcJ = problem_.variables().capillaryPressure(
                            globalIdxJ);

                    // calculate the saturation gradient
                    Dune::FieldVector<Scalar, dimWorld> pcGradient =
                            unitOuterNormal;
                    pcGradient *= (pcI - pcJ) / dist;

                    // get the diffusive part
                    Scalar diffPart = diffusivePart()(*eIt, isIndex, satWI,
                            satWJ, pcGradient) * unitOuterNormal * faceArea;


                    Scalar convPart = convectivePart()(*eIt, isIndex, satWI,
                            satWJ) * unitOuterNormal * faceArea;

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
                        diffPart *= -1;
                        convPart *= -1;
                        break;
                    }
                    }
                    factor -= diffPart;
                    factor += convPart;

                    //add cflFlux for time-stepping
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, 10 * diffPart, *isIt);
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, 10 * convPart, *isIt);

                    break;
                }
                case vw:
                {
                    if (compressibility_)
                    {
                        factor /= densityWI;
                        factorSecondPhase /= densityNWI;
                    }

                    //add cflFlux for time-stepping
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, factor, *isIt, wPhaseIdx);
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, factorSecondPhase,
                            *isIt, nPhaseIdx);

                    break;
                }
                case vn:
                {
                    if (compressibility_)
                    {
                        factor /= densityNWI;
                        factorSecondPhase /= densityWI;
                    }

                    //add cflFlux for time-stepping
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, factor, *isIt, nPhaseIdx);
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                            viscosityWI, viscosityNWI, factorSecondPhase,
                            *isIt, wPhaseIdx);

                    break;
                }
                }
            }//end intersection with neighbor element

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector<Scalar, dimWorld> distVec = globalPosFace
                        - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                //get boundary type
                BoundaryConditions::Flags bcTypeSat = problem_.bctypeSat(
                        globalPosFace, *isIt);

                if (bcTypeSat == BoundaryConditions::dirichlet)
                {
                    Scalar satBound = problem_.dirichletSat(globalPosFace,
                            *isIt);

                    //get velocity*normalvector*facearea/(volume*porosity)
                    factor
                            = (problem_.variables().velocity()[globalIdxI][isIndex]
                                    * unitOuterNormal) * (faceArea);
                    factorSecondPhase
                            = (problem_.variables().velocitySecondPhase()[globalIdxI][isIndex]
                                    * unitOuterNormal) * (faceArea);

                    Scalar pressBound =
                            problem_.variables().pressure()[globalIdxI];
                    Scalar temperature = problem_.temperature(globalPosFace,
                            *eIt);

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
                        satWI = 1
                                - problem_.variables().saturation()[globalIdxI];
                        satWBound = 1 - satBound;
                        break;
                    }
                    }

                    Scalar pcBound = MaterialLaw::pC(
                            problem_.spatialParameters().materialLawParams(
                                    globalPos, *eIt), satBound);

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
                    Scalar potentialW = problem_.variables().potentialWetting(
                            globalIdxI, isIndex);
                    Scalar potentialNW =
                            problem_.variables().potentialNonwetting(
                                    globalIdxI, isIndex);

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
                            lambdaW
                                    = MaterialLaw::krw(
                                            problem_.spatialParameters().materialLawParams(
                                                    globalPos, *eIt), satWBound)
                                            / FluidSystem::phaseViscosity(
                                                    wPhaseIdx, temperature,
                                                    pressW, fluidState);
                        }
                        else
                        {
                            lambdaW
                                    = MaterialLaw::krw(
                                            problem_.spatialParameters().materialLawParams(
                                                    globalPos, *eIt), satWBound)
                                            / viscosityWI;
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
                            lambdaNW
                                    = MaterialLaw::krn(
                                            problem_.spatialParameters().materialLawParams(
                                                    globalPos, *eIt), satWBound)
                                            / FluidSystem::phaseViscosity(
                                                    nPhaseIdx, temperature,
                                                    pressNW, fluidState);
                        }
                        else
                        {
                            lambdaNW
                                    = MaterialLaw::krn(
                                            problem_.spatialParameters().materialLawParams(
                                                    globalPos, *eIt), satWBound)
                                            / viscosityNWI;
                        }
                    }
                    //                    std::cout<<lambdaW<<" "<<lambdaNW<<std::endl;

                    switch (velocityType_)
                    {
                    case vt:
                    {
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt);

                        Scalar pcI = problem_.variables().capillaryPressure(
                                globalIdxI);

                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar, dimWorld> pcGradient =
                                unitOuterNormal;
                        pcGradient *= (pcI - pcBound) / dist;

                        // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
                        Scalar diffPart = diffusivePart()(*eIt, isIndex, satWI,
                                satWBound, pcGradient) * unitOuterNormal
                                * faceArea;

                        Scalar convPart = convectivePart()(*eIt, isIndex,
                                satWI, satWBound) * unitOuterNormal * faceArea;

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
                            diffPart *= -1;
                            convPart *= -1;
                            break;
                        }
                        }
                        //vt*fw
                        factor -= diffPart;
                        factor += convPart;

                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, 10 * diffPart, *isIt);
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, 10 * convPart, *isIt);

                        break;
                    }

                    case vw:
                    {
                        if (compressibility_)
                        {
                            factor /= densityWI;
                            factorSecondPhase /= densityNWI;
                        }

                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt,
                                wPhaseIdx);
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factorSecondPhase,
                                *isIt, nPhaseIdx);

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

                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt,
                                nPhaseIdx);
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factorSecondPhase,
                                *isIt, wPhaseIdx);

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

                    switch (saturationType_)
                    {
                    case Sw:
                    {
                        factor = problem_.neumann(globalPosFace, *isIt)[wPhaseIdx];
                        factor /= densityWI * faceArea;
                        break;
                    }
                    case Sn:
                    {
                        factor = problem_.neumann(globalPosFace, *isIt)[nPhaseIdx];
                        factor /= densityNWI * faceArea;
                        break;
                    }
                    }

                    switch (velocityType_)
                    {
                    case vt:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt);
                    }
                    case vw:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt,
                                wPhaseIdx);
                        break;
                    }
                    case vn:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt,
                                nPhaseIdx);
                        break;
                    }
                    }

                }//end neumann boundary
                if (bcTypeSat == BoundaryConditions::outflow)
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

                    //get velocity*normalvector*facearea/(volume*porosity)
                    factor
                            = (problem_.variables().velocity()[globalIdxI][isIndex]
                                    * unitOuterNormal) * faceArea;

                    if (velocityType_ == vt)
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
                    }

                    switch (velocityType_)
                    {
                    case vt:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt);
                    }
                    case vw:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt,
                                wPhaseIdx);
                        break;
                    }
                    case vn:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW,
                                viscosityWI, viscosityNWI, factor, *isIt,
                                nPhaseIdx);
                        break;
                    }
                    }

                }//end outflow boundary
            }//end boundary
            // add to update vector
            updateVec[globalIdxI] -= factor / (volume * porosity);//-:v>0, if flow leaves the cell
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
            source = problem_.source(globalPos, *eIt)[wPhaseIdx] / densityWI
                    + problem_.source(globalPos, *eIt)[nPhaseIdx] / densityNWI;
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

            switch (velocityType_)
            {
            case vw:
            {
                //add cflFlux for time-stepping
                evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityWI,
                        viscosityNWI, source * volume, *eIt, wPhaseIdx);
                break;
            }
            case vn:
            {
                //add cflFlux for time-stepping
                evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityWI,
                        viscosityNWI, source * volume, *eIt, nPhaseIdx);
                break;
            }
            case vt:
            {
                //add cflFlux for time-stepping
                evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityWI,
                        viscosityNWI, source * volume, *eIt);
                break;
            }
            }
        }

        //calculate time step
        dt = std::min(dt, evalCflFluxFunction().getCflFluxFunction(globalPos,
                *eIt) * (porosity * volume));

        problem_.variables().volumecorrection(globalIdxI)
                = updateVec[globalIdxI];
    } // end grid traversal

    return 0;
}

template<class TypeTag>
void FVSaturation2P<TypeTag>::initialize()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt
            != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        // initialize cell concentration
        problem_.variables().saturation()[problem_.variables().index(*eIt)]
                = problem_.initSat(globalPos, *eIt);
    }
    return;
}

template<class TypeTag>
void FVSaturation2P<TypeTag>::updateMaterialLaws(
        RepresentationType& saturation = *(new RepresentationType(0)),
        bool iterate = false)
{
    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt
            != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

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
        Scalar referencePressure = problem_.referencePressure(globalPos, *eIt);

        fluidState.update(satW, referencePressure, referencePressure,
                temperature);

        problem_.variables().densityWetting(globalIdx)
                = FluidSystem::phaseDensity(wPhaseIdx, temperature,
                        referencePressure, fluidState);
        problem_.variables().densityNonwetting(globalIdx)
                = FluidSystem::phaseDensity(nPhaseIdx, temperature,
                        referencePressure, fluidState);
        problem_.variables().viscosityWetting(globalIdx)
                = FluidSystem::phaseViscosity(wPhaseIdx, temperature,
                        referencePressure, fluidState);
        problem_.variables().viscosityNonwetting(globalIdx)
                = FluidSystem::phaseViscosity(nPhaseIdx, temperature,
                        referencePressure, fluidState);

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx)
                = MaterialLaw::krw(
                        problem_.spatialParameters().materialLawParams(
                                globalPos, *eIt), satW)
                        / problem_.variables().viscosityWetting(globalIdx);
        problem_.variables().mobilityNonwetting(globalIdx)
                = MaterialLaw::krn(
                        problem_.spatialParameters().materialLawParams(
                                globalPos, *eIt), satW)
                        / problem_.variables().viscosityNonwetting(globalIdx);
        problem_.variables().capillaryPressure(globalIdx)
                = MaterialLaw::pC(
                        problem_.spatialParameters().materialLawParams(
                                globalPos, *eIt), satW);

        problem_.variables().fracFlowFuncWetting(globalIdx)
                = problem_.variables().mobilityWetting(globalIdx)
                        / (problem_.variables().mobilityWetting(globalIdx)
                                + problem_.variables().mobilityNonwetting(
                                        globalIdx));
        problem_.variables().fracFlowFuncNonwetting(globalIdx)
                = problem_.variables().mobilityNonwetting(globalIdx)
                        / (problem_.variables().mobilityWetting(globalIdx)
                                + problem_.variables().mobilityNonwetting(
                                        globalIdx));

        problem_.spatialParameters().update(satW, *eIt);
    }
    return;
}

}
#endif
