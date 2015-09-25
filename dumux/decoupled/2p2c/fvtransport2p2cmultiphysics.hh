// $Id:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff, Benjamin Faigle                     *
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
#ifndef DUMUX_FVTRANSPORT2P2CMULTIPHYSICS_HH
#define DUMUX_FVTRANSPORT2P2CMULTIPHYSICS_HH

#include <dumux/decoupled/2p2c/2p2cproperties.hh>
#include <dumux/common/math.hh>

/**
 * @file
 * @brief  Finite Volume discretization of the component transport equation
 * @author Markus Wolff, Jochen Fritz, Benjamin Faigle
 */

namespace Dumux
{
//! Miscible Transport step in a Finite Volume discretization
/*!
 * \ingroup multiphysics
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *
 * The model domain is automatically divided
 * in a single-phase and a two-phase domain. The full 2p2c model is only evaluated within the
 * two-phase subdomain, whereas a single-phase transport model is computed in the rest of the
 * domain.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2CMultiPhysics
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TransportSolutionType)) TransportSolutionType;

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
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes);

    //! Set the initial values before the first pressure equation
    /*!
     * This method is called before first pressure equation is solved from Dumux::IMPET.
     */
    void initialize()
    {};

    //! \brief Write data files
     /*  \param writer applied VTK-writer */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        return;
    }


    //! Constructs a FVTransport2P2CMultiPhysics object
    /**
     * \param problem a problem class object
     */

    FVTransport2P2CMultiPhysics(Problem& problem) :
        problem_(problem), switchNormals(false)
    {
        const int velocityType = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation));
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility)) && velocityType == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        const int saturationType = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation));
        if (saturationType != Sw && saturationType != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (pressureType != pw && pressureType != pn && pressureType != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (velocityType != vw && velocityType != vn && velocityType != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }
    }

    ~FVTransport2P2CMultiPhysics()
    {
    }

private:
    Problem& problem_;

    static const int pressureType = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    bool switchNormals;
};

//! \brief Calculate the update vector and determine timestep size
/*!
 *  This method calculates the update vector \f$ u \f$ of the discretized equation
 *  \f[
       C^{\kappa , new} = C^{\kappa , old} + u,
 *  \f]
 *  where \f$ u = \sum_{element faces} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{element face} \f$,
 *  \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{element face} \f$ is the face area.
 *
 *  In addition to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 *  This method = old concentrationUpdate()
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param[out] dt Time step size \f$\mathrm{[s]}\f$
 *  \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2CMultiPhysics<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false)
{
    // initialize dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    // some phase properties
    Dune::FieldVector<Scalar, dimWorld> gravity_ = problem_.gravity();

    // compute update vector
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        if(impet || (problem_.variables().subdomain(globalIdxI)==2))   // estimate only necessary in subdomain
        {
            // cell geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get position
            GlobalPosition globalPos = eIt->geometry().center();

            // cell volume, assume linear map here
            Scalar volume = eIt->geometry().volume();

            // get values of cell I
            Scalar pressI = problem_.variables().pressure()[globalIdxI];
    //        Scalar pcI = problem_.variables().capillaryPressure(globalIdxI);
            Dune::FieldMatrix<Scalar,dim,dim> K_I(problem_.spatialParameters().intrinsicPermeability(globalPos, *eIt));

            Scalar SwmobI = std::max((problem_.variables().saturation(globalIdxI)
                                    - problem_.spatialParameters().materialLawParams(globalPos, *eIt).Swr())
                                    , 1e-2);
            Scalar SnmobI = std::max((1. - problem_.variables().saturation(globalIdxI)
                                        - problem_.spatialParameters().materialLawParams(globalPos, *eIt).Snr())
                                    , 1e-2);

            double Xw1_I = problem_.variables().wet_X1(globalIdxI);
            double Xn1_I = problem_.variables().nonwet_X1(globalIdxI);

            Scalar densityWI = problem_.variables().densityWetting(globalIdxI);
            Scalar densityNWI = problem_.variables().densityNonwetting(globalIdxI);

            // some variables for time step calculation
            double sumfactorin = 0;
            double sumfactorout = 0;

            // run through all intersections with neighbors and boundary
            IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                // get geometry type of face
                Dune::GeometryType faceGT = isIt->geometryInInside().type();

                // center in face's reference element
                typedef Dune::GenericReferenceElements<Scalar, dim - 1> FaceReferenceElements;
                const Dune::FieldVector<Scalar, dim - 1>& faceLocal =
                    FaceReferenceElements::general(faceGT).position(0, 0);

                // center of face inside volume reference element
                const LocalPosition localPosFace(0);

                Dune::FieldVector<Scalar, dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);
                if (switchNormals)
                    unitOuterNormal *= -1.0;

                Scalar faceArea = isIt->geometry().volume();

                // create vector for timestep and for update
                Dune::FieldVector<Scalar, 2> factor (0.);
                Dune::FieldVector<Scalar, 2> updFactor (0.);

                // handle interior face
                if (isIt->neighbor())
                {
                    // access neighbor
                    ElementPointer neighborPointer = isIt->outside();
                    int globalIdxJ = problem_.variables().index(*neighborPointer);

                    // compute factor in neighbor
                    Dune::GeometryType neighborGT = neighborPointer->geometry().type();

                    // cell center in global coordinates
                    const GlobalPosition& globalPos = eIt->geometry().center();

                    // neighbor cell center in global coordinates
                    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar, dimWorld> distVec = globalPosNeighbor - globalPos;
                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                    unitDistVec /= dist;

                    // get saturation and concentration value at neighbor cell center
                    double Xw1_J = problem_.variables().wet_X1(globalIdxJ);
                    double Xn1_J = problem_.variables().nonwet_X1(globalIdxJ);

                    // phase densities in neighbor
                    Scalar densityWJ = problem_.variables().densityWetting(globalIdxJ);
                    Scalar densityNWJ = problem_.variables().densityNonwetting(globalIdxJ);

                    // average phase densities with central weighting
                    double densityW_mean = (densityWI + densityWJ) * 0.5;
                    double densityNW_mean = (densityNWI + densityNWJ) * 0.5;

                    double pressJ = problem_.variables().pressure()[globalIdxJ];

                    // compute mean permeability
                    Dune::FieldMatrix<Scalar,dim,dim> meanK_(0.);
                    Dumux::harmonicMeanMatrix(meanK_,
                            K_I,
                            problem_.spatialParameters().intrinsicPermeability(globalPosNeighbor, *neighborPointer));
                    Dune::FieldVector<Scalar,dim> K(0);
                    meanK_.umv(unitDistVec,K);

                    // velocities
                    double potentialW = (K * unitOuterNormal) * (pressI - pressJ) / (dist);
                    double potentialNW = potentialW + (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityNW_mean;
                    potentialW += (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityW_mean;
                    //                  lambda = 2 * lambdaI * lambdaJ / (lambdaI + lambdaJ);

                    double lambdaW, lambdaNW;

                    if (potentialW >= 0.)
                        lambdaW = problem_.variables().mobilityWetting(globalIdxI);
                    else
                        lambdaW = problem_.variables().mobilityWetting(globalIdxJ);

                    if (potentialNW >= 0.)
                        lambdaNW = problem_.variables().mobilityNonwetting(globalIdxI);
                    else
                        lambdaNW = problem_.variables().mobilityNonwetting(globalIdxJ);

                    // calculate and standardized velocity
                    double velocityJIw = std::max(-lambdaW * potentialW * faceArea / volume, 0.0);
                    double velocityIJw = std::max( lambdaW * potentialW * faceArea / volume, 0.0);
                    double velocityJIn = std::max(-lambdaNW * potentialNW * faceArea / volume, 0.0);
                    double velocityIJn = std::max( lambdaNW * potentialNW * faceArea / volume, 0.0);

                    // for timestep control
                    factor[0] = velocityJIw + velocityJIn;

                    double foutw = velocityIJw/SwmobI;
                    double foutn = velocityIJn/SnmobI;
                    if (std::isnan(foutw) || std::isinf(foutw) || foutw < 0) foutw = 0;
                    if (std::isnan(foutn) || std::isinf(foutn) || foutn < 0) foutn = 0;
                    factor[1] = foutw + foutn;

                    updFactor[wCompIdx] =
                        velocityJIw * Xw1_J * densityWJ
                        - velocityIJw * Xw1_I * densityWI
                        + velocityJIn * Xn1_J * densityNWJ
                        - velocityIJn * Xn1_I * densityNWI;
                    updFactor[nCompIdx] =
                        velocityJIw * (1. - Xw1_J) * densityWJ
                        - velocityIJw * (1. - Xw1_I) * densityWI
                        + velocityJIn * (1. - Xn1_J) * densityNWJ
                        - velocityIJn * (1. - Xn1_I) * densityNWI;
                }

                /******************************************
                 *     Boundary Face
                 ******************************************/
                if (isIt->boundary())
                {
                    // center of face in global coordinates
                    GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                    // cell center in global coordinates
                    GlobalPosition globalPos = eIt->geometry().center();

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar, dimWorld> distVec = globalPosFace - globalPos;
                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                    unitDistVec /= dist;

                    //instantiate a fluid state
                    FluidState fluidState;

                    //get boundary type
                    BoundaryConditions::Flags bcTypeTransport_ = problem_.bcTypeTransport(globalPosFace, *isIt);

                    if (bcTypeTransport_ == BoundaryConditions::dirichlet)
                    {
                        //get dirichlet pressure boundary condition
                        Scalar pressBound = problem_.dirichletPress(globalPosFace, *isIt);

                        // read boundary values
                        BoundaryConditions2p2c::Flags bctype = problem_.bcFormulation(globalPosFace, *isIt);
                        if (bctype == BoundaryConditions2p2c::saturation)
                        {
                            Scalar satBound = problem_.dirichletTransport(globalPosFace, *isIt);
                            fluidState.satFlash(satBound, pressBound,
                                                problem_.spatialParameters().porosity(globalPos, *eIt),
                                                problem_.temperature(globalPosFace, *eIt));

                        }
                        if (bctype == BoundaryConditions2p2c::concentration)
                        {
                            Scalar Z1Bound = problem_.dirichletTransport(globalPosFace, *isIt);
                            fluidState.update(Z1Bound, pressBound,
                                                problem_.spatialParameters().porosity(globalPos, *eIt),
                                                problem_.temperature(globalPosFace, *eIt));
                        }
                        // determine fluid properties at the boundary
                        Scalar Xw1Bound = fluidState.massFrac(wPhaseIdx, wCompIdx);
                        Scalar Xn1Bound = fluidState.massFrac(nPhaseIdx, wCompIdx);
                        Scalar densityWBound = FluidSystem::phaseDensity(wPhaseIdx,
                                                                    problem_.temperature(globalPosFace, *eIt),
                                                                    pressBound, fluidState);
                        Scalar densityNWBound = FluidSystem::phaseDensity(nPhaseIdx,
                                                                    problem_.temperature(globalPosFace, *eIt),
                                                                    pressBound, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx,
                                                                            problem_.temperature(globalPosFace, *eIt),
                                                                            pressBound, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx,
                                                                            problem_.temperature(globalPosFace, *eIt),
                                                                            pressBound, fluidState);

                        // average
                        double densityW_mean = (densityWI + densityWBound) / 2;
                        double densityNW_mean = (densityNWI + densityNWBound) / 2;

                        // prepare K
                        Dune::FieldVector<Scalar,dim> K(0);
                        K_I.umv(unitDistVec,K);

                        // velocities
                        double potentialW = (K * unitOuterNormal) * (pressI - pressBound) / (dist);
                        double potentialNW = potentialW + (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityNW_mean;
                        potentialW += (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityW_mean;

                        double lambdaW, lambdaN;

                        if (potentialW >= 0.)
                            lambdaW = problem_.variables().mobilityWetting(globalIdxI);
                        else
                            {
                            if(GET_PROP_VALUE(TypeTag, PTAG(BoundaryMobility))==Indices::satDependent)
                                lambdaW = fluidState.saturation(wPhaseIdx) / viscosityWBound;
                            else
                                lambdaW = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                                        / viscosityWBound;
                            }
                        if (potentialNW >= 0.)
                            lambdaN = problem_.variables().mobilityNonwetting(globalIdxI);
                        else
                            {
                            if(GET_PROP_VALUE(TypeTag, PTAG(BoundaryMobility))==Indices::satDependent)
                                lambdaN = fluidState.saturation(nPhaseIdx) / viscosityNWBound;
                            else
                                lambdaN = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(globalPos, *eIt), fluidState.saturation(wPhaseIdx))
                                        / viscosityNWBound;
                            }
                        // standardized velocity
                        double velocityJIw = std::max(-lambdaW * potentialW * faceArea / volume, 0.0);
                        double velocityIJw = std::max( lambdaW * potentialW * faceArea / volume, 0.0);
                        double velocityJIn = std::max(-lambdaN * potentialNW * faceArea / volume, 0.0);
                        double velocityIJn = std::max( lambdaN * potentialNW* faceArea / volume, 0.0);

                        // for timestep control
                        factor[0] = velocityJIw + velocityJIn;

                        double foutw = velocityIJw/SwmobI;
                        double foutn = velocityIJn/SnmobI;
                        if (std::isnan(foutw) || std::isinf(foutw) || foutw < 0) foutw = 0;
                        if (std::isnan(foutn) || std::isinf(foutn) || foutn < 0) foutn = 0;
                        factor[1] = foutw + foutn;

                        updFactor[wCompIdx] =
                            + velocityJIw * Xw1Bound * densityWBound
                            - velocityIJw * Xw1_I * densityWI
                            + velocityJIn * Xn1Bound * densityNWBound
                            - velocityIJn * Xn1_I * densityNWI ;
                        updFactor[nCompIdx] =
                            velocityJIw * (1. - Xw1Bound) * densityWBound
                            - velocityIJw * (1. - Xw1_I) * densityWI
                            + velocityJIn * (1. - Xn1Bound) * densityNWBound
                            - velocityIJn * (1. - Xn1_I) * densityNWI ;
                    }//end dirichlet boundary

                    if (bcTypeTransport_ == BoundaryConditions::neumann)
                    {
                        // Convention: outflow => positive sign : has to be subtracted from update vec
                        Dune::FieldVector<Scalar,2> J = problem_.neumann(globalPosFace, *isIt);
                        updFactor[wCompIdx] = - J[0] * faceArea / volume;
                        updFactor[nCompIdx] = - J[1] * faceArea / volume;

                        // for timestep control
                        #define cflIgnoresNeumann
                        #ifdef cflIgnoresNeumann
                        factor[0] = 0;
                        factor[1] = 0;
                        #else
                        double inflow = updFactor[wCompIdx] / densityW + updFactor[nCompIdx] / densityNW;
                        if (inflow>0)
                            {
                            factor[0] = updFactor[wCompIdx] / densityW + updFactor[nCompIdx] / densityNW;    // =factor in
                            factor[1] = -(updFactor[wCompIdx] / densityW /SwmobI + updFactor[nCompIdx] / densityNW / SnmobI);    // =factor out
                            }
                        else
                        {
                            factor[0] = -(updFactor[wCompIdx] / densityW + updFactor[nCompIdx] / densityNW);    // =factor in
                            factor[1] = updFactor[wCompIdx] / densityW /SwmobI + updFactor[nCompIdx] / densityNW / SnmobI;    // =factor out
                        }
                        #endif
                    }//end neumann boundary
                }//end boundary
                // add to update vector
                updateVec[wCompIdx][globalIdxI] += updFactor[wCompIdx];
                updateVec[nCompIdx][globalIdxI] += updFactor[nCompIdx];

                // for time step calculation
                sumfactorin += factor[0];
                sumfactorout += factor[1];

            }// end all intersections

            /************************************
             *     Handle source term
             ***********************************/
            Dune::FieldVector<double,2> q = problem_.source(globalPos, *eIt);
            updateVec[wCompIdx][globalIdxI] += q[wCompIdx];
            updateVec[nCompIdx][globalIdxI] += q[nCompIdx];

            // account for porosity
            sumfactorin = std::max(sumfactorin,sumfactorout)
                            / problem_.spatialParameters().porosity(globalPos, *eIt);

            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= globalIdxI;
            }
        }
    } // end grid traversal
    if(impet)
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "<<dt * GET_PROP_VALUE(TypeTag, PTAG(CFLFactor))<< std::endl;

    return;
}

}
#endif
