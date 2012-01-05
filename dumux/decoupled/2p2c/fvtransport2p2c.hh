// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef DUMUX_FVTRANSPORT2P2C_HH
#define DUMUX_FVTRANSPORT2P2C_HH

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
/*! \ingroup multiphase
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2C
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

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
        Sn = Indices::saturationNW
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

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, 2> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    //! Acess function for the current problem
    Problem& problem()
    {return problem_;};
public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes);

    //! Set the initial values before the first pressure equation
    /*!
     * This method is called before first pressure equation is solved from Dumux::IMPET.
     */
    void initialize()
    {};

    void evalBoundary(GlobalPosition globalPosFace, IntersectionIterator& isIt,
                         ElementIterator& eIt, FluidState &BCfluidState, PhaseVector &pressure);

    //! \brief Write data files
     /*  \param writer applied VTK-writer */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        return;
    }


    //! Constructs a FVTransport2P2C object
    /*!
     * Currently, the miscible transport scheme can not be applied with a global pressure / total velocity
     * formulation.
     *
     * \param problem a problem class object
     */
    FVTransport2P2C(Problem& problem) :
        problem_(problem), switchNormals(false)
    {  }

    virtual ~FVTransport2P2C()
    {
    }

protected:
    Problem& problem_;

    static const int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
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
void FVTransport2P2C<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false)
{
    // initialize dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    // some phase properties
    const GlobalPosition& gravity_ = problem().gravity();

    // compute update vector
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
#if HAVE_MPI
        if (eIt->partitionType() != Dune::InteriorEntity)
        {
            continue;
        }
#endif

        // cell index
        int globalIdxI = problem().variables().index(*eIt);

        // get position
        GlobalPosition globalPos = eIt->geometry().center();

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        // get values of cell I
        Scalar pressI = problem().variables().pressure()[globalIdxI];
        Scalar pcI = problem().variables().capillaryPressure(globalIdxI);
        Dune::FieldMatrix<Scalar,dim,dim> K_I(problem().spatialParameters().intrinsicPermeability(globalPos, *eIt));

        Scalar SwmobI = std::max((problem().variables().saturation(globalIdxI)
                                - problem().spatialParameters().materialLawParams(globalPos, *eIt).Swr())
                                , 1e-2);
        Scalar SnmobI = std::max((1. - problem().variables().saturation(globalIdxI)
                                    - problem().spatialParameters().materialLawParams(globalPos, *eIt).Snr())
                                , 1e-2);


        double Xw1_I = problem().variables().wet_X1(globalIdxI);
        double Xn1_I = problem().variables().nonwet_X1(globalIdxI);

        Scalar densityWI (0.), densityNWI(0.);
        if (GET_PROP_VALUE(TypeTag, NumDensityTransport))
        {
            densityWI= problem().variables().numericalDensity(globalIdxI, wPhaseIdx);
            densityNWI = problem().variables().numericalDensity(globalIdxI, nPhaseIdx);

        }
        else
        {
            densityWI= problem().variables().densityWetting(globalIdxI);
            densityNWI = problem().variables().densityNonwetting(globalIdxI);
        }
        // some variables for time step calculation
        double sumfactorin = 0;
        double sumfactorout = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            GlobalPosition unitOuterNormal = isIt->centerUnitOuterNormal();
            if (switchNormals)
                unitOuterNormal *= -1.0;

            Scalar faceArea = isIt->geometry().volume();

            // create vector for timestep and for update
            Dune::FieldVector<Scalar, 2> factor (0.);
            Dune::FieldVector<Scalar, 2> updFactor (0.);

            Scalar potentialW(0.), potentialNW(0.);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem().variables().index(*neighborPointer);

                // cell center in global coordinates
                const GlobalPosition& globalPos = eIt->geometry().center();

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosNeighbor - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                // get saturation and concentration value at neighbor cell center
                double Xw1_J = problem().variables().wet_X1(globalIdxJ);
                double Xn1_J = problem().variables().nonwet_X1(globalIdxJ);

                // phase densities in neighbor
                Scalar densityWJ (0.), densityNWJ(0.);
                if (GET_PROP_VALUE(TypeTag, NumDensityTransport))
                {
                    densityWJ = problem().variables().numericalDensity(globalIdxJ, wPhaseIdx);
                    densityNWJ = problem().variables().numericalDensity(globalIdxJ, nPhaseIdx);
                }
                else
                {
                    densityWJ = problem().variables().densityWetting(globalIdxJ);
                    densityNWJ = problem().variables().densityNonwetting(globalIdxJ);
                }

                // average phase densities with central weighting
                double densityW_mean = (densityWI + densityWJ) * 0.5;
                double densityNW_mean = (densityNWI + densityNWJ) * 0.5;

                double pressJ = problem().variables().pressure()[globalIdxJ];
                Scalar pcJ = problem().variables().capillaryPressure(globalIdxJ);

                // compute mean permeability
                Dune::FieldMatrix<Scalar,dim,dim> meanK_(0.);
                Dumux::harmonicMeanMatrix(meanK_,
                        K_I,
                        problem().spatialParameters().intrinsicPermeability(globalPosNeighbor, *neighborPointer));
                Dune::FieldVector<Scalar,dim> K(0);
                meanK_.umv(unitDistVec,K);

                // determine potentials for upwind
                switch (pressureType)
                {
                case pw:
                {
                    potentialW = (K * unitOuterNormal) * (pressI - pressJ) / (dist);
                    potentialNW = (K * unitOuterNormal) * (pressI - pressJ + pcI - pcJ) / (dist);
                    break;
                }
                case pn:
                {
                    potentialW = (K * unitOuterNormal) * (pressI - pressJ - pcI + pcJ) / (dist);
                    potentialNW = (K * unitOuterNormal) * (pressI - pressJ) / (dist);
                    break;
                }
                }
                // add gravity term
                potentialNW +=  (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityNW_mean;
                potentialW +=  (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityW_mean;

                // upwind mobility
                double lambdaW, lambdaNW;
                if (potentialW >= 0.)
                    lambdaW = problem().variables().mobilityWetting(globalIdxI);
                else
                    lambdaW = problem().variables().mobilityWetting(globalIdxJ);

                if (potentialNW >= 0.)
                    lambdaNW = problem().variables().mobilityNonwetting(globalIdxI);
                else
                    lambdaNW = problem().variables().mobilityNonwetting(globalIdxJ);

                // calculate and standardized velocity
                double velocityJIw = std::max((-lambdaW * potentialW) * faceArea / volume, 0.0);
                double velocityIJw = std::max(( lambdaW * potentialW) * faceArea / volume, 0.0);
                double velocityJIn = std::max((-lambdaNW * potentialNW) * faceArea / volume, 0.0);
                double velocityIJn = std::max(( lambdaNW * potentialNW) * faceArea / volume, 0.0);

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
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                // cell center in global coordinates
                const GlobalPosition& globalPos = eIt->geometry().center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosFace - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                //instantiate a fluid state
                FluidState BCfluidState;

                //get boundary type
                BoundaryTypes bcTypes;
                problem().boundaryTypes(bcTypes, *isIt);

                /**********         Dirichlet Boundary        *************/
                if (bcTypes.isDirichlet(Indices::contiWEqIdx)) // if contiWEq is Dirichlet, so is contiNEq
                {
                    //get dirichlet pressure boundary condition
                    PhaseVector pressBound(0.);
                    Scalar pcBound (0.);

                    // read boundary values
                    this->evalBoundary(globalPosFace,
                                    isIt,
                                    eIt,
                                    BCfluidState,
                                    pressBound);

                    // determine fluid properties at the boundary
                    Scalar Xw1Bound = BCfluidState.massFraction(wPhaseIdx, wCompIdx);
                    Scalar Xn1Bound = BCfluidState.massFraction(nPhaseIdx, wCompIdx);
                    Scalar densityWBound = BCfluidState.density(wPhaseIdx);
                    Scalar densityNWBound = BCfluidState.density(nPhaseIdx);
                    Scalar viscosityWBound = FluidSystem::viscosity(BCfluidState, wPhaseIdx);
                    Scalar viscosityNWBound = FluidSystem::viscosity(BCfluidState, nPhaseIdx);
                    if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
                        pcBound = BCfluidState.capillaryPressure();
                    // average
                    double densityW_mean = (densityWI + densityWBound) / 2;
                    double densityNW_mean = (densityNWI + densityNWBound) / 2;

                    // prepare K
                    Dune::FieldVector<Scalar,dim> K(0);
                    K_I.umv(unitDistVec,K);

                    //calculate potential gradient
                    switch (pressureType)
                    {
                        case pw:
                        {
                            potentialW = (K * unitOuterNormal) *
                                    (pressI - pressBound[wPhaseIdx]) / dist;
                            potentialNW = (K * unitOuterNormal) *
                                    (pressI + pcI - pressBound[wPhaseIdx] - pcBound)
                                    / dist;
                            break;
                        }
                        case pn:
                        {
                            potentialW = (K * unitOuterNormal) *
                                    (pressI - pcI - pressBound[nPhaseIdx] + pcBound)
                                    / dist;
                            potentialNW = (K * unitOuterNormal) *
                                    (pressI - pressBound[nPhaseIdx]) / dist;
                            break;
                        }
                    }
                    potentialW += (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityW_mean;;
                    potentialNW += (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityNW_mean;;

                    // do upwinding for lambdas
                    double lambdaW, lambdaNW;
                    if (potentialW >= 0.)
                        lambdaW = problem().variables().mobilityWetting(globalIdxI);
                    else
                        {
                        if(GET_PROP_VALUE(TypeTag, BoundaryMobility)==Indices::satDependent)
                            lambdaW = BCfluidState.saturation(wPhaseIdx) / viscosityWBound;
                        else
                            lambdaW = MaterialLaw::krw(
                                    problem().spatialParameters().materialLawParams(globalPos, *eIt), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityWBound;
                        }
                    if (potentialNW >= 0.)
                        lambdaNW = problem().variables().mobilityNonwetting(globalIdxI);
                    else
                        {
                        if(GET_PROP_VALUE(TypeTag, BoundaryMobility)==Indices::satDependent)
                            lambdaNW = BCfluidState.saturation(nPhaseIdx) / viscosityNWBound;
                        else
                            lambdaNW = MaterialLaw::krn(
                                    problem().spatialParameters().materialLawParams(globalPos, *eIt), BCfluidState.saturation(wPhaseIdx))
                                    / viscosityNWBound;
                        }
                    // calculate and standardized velocity
                    double velocityJIw = std::max((-lambdaW * potentialW) * faceArea / volume, 0.0);
                    double velocityIJw = std::max(( lambdaW * potentialW) * faceArea / volume, 0.0);
                    double velocityJIn = std::max((-lambdaNW * potentialNW) * faceArea / volume, 0.0);
                    double velocityIJn = std::max(( lambdaNW * potentialNW) * faceArea / volume, 0.0);

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
                else if (bcTypes.isNeumann(Indices::contiWEqIdx))
                {
                    // Convention: outflow => positive sign : has to be subtracted from update vec
                    PrimaryVariables J(NAN);
                    problem().neumann(J, *isIt);
                    updFactor[wCompIdx] = - J[Indices::contiWEqIdx] * faceArea / volume;
                    updFactor[nCompIdx] = - J[Indices::contiNEqIdx] * faceArea / volume;

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
            // correct update Factor by volume error
            #ifdef errorInTransport
            updFactor[wCompIdx] *= (1+problem().variables().volErr(globalIdxI));
            updFactor[nCompIdx] *= (1+problem().variables().volErr(globalIdxI));
            #endif
            // add to update vector
            updateVec[wCompIdx][globalIdxI] += updFactor[wCompIdx];
            updateVec[nCompIdx][globalIdxI] += updFactor[nCompIdx];

            // for time step calculation
            sumfactorin += factor[0];
            sumfactorout += factor[1];

        }// end all intersections

        /***********     Handle source term     ***************/
        PrimaryVariables q(NAN);
        problem().source(q, *eIt);
        updateVec[wCompIdx][globalIdxI] += q[Indices::contiWEqIdx];
        updateVec[nCompIdx][globalIdxI] += q[Indices::contiNEqIdx];

        // account for porosity
        sumfactorin = std::max(sumfactorin,sumfactorout)
                        / problem().spatialParameters().porosity(globalPos, *eIt);

        if ( 1./sumfactorin < dt)
        {
            dt = 1./sumfactorin;
            restrictingCell= globalIdxI;
        }
    } // end grid traversal
    if(impet)
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "<<dt * GET_PARAM(TypeTag, Scalar, CFLFactor)<< std::endl;
    return;
}

//! evaluate the boundary conditions
/*!
 *  As the transport primary variable in this formulation is the total component concentration, \f$ C^{\kappa} \f$
 *  it seems natural that the boundary values are also total concentrations. However, as for the initial conditions,
 *  it is possible to define boundaries by means of a saturation. This choice determines which version of flash
 *  calculation is necessary to get to the composition at the boundary.
 *
 *  \param globalPosFace Face of the current boundary
 *  \param isIt Iterator pointing on the current intersection
 *  \param eIt Iterator pointing on the current element
 *  \param BCfluidState FluidState object that is used for the boundary
 *  \param pressBound Boundary values of phase pressures
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::evalBoundary(GlobalPosition globalPosFace,
                                            IntersectionIterator& isIt,
                                            ElementIterator& eIt,
                                            FluidState &BCfluidState,
                                            PhaseVector &pressBound)
{
    // read boundary values
    PrimaryVariables primaryVariablesOnBoundary(0.);
    problem().dirichlet(primaryVariablesOnBoundary, *isIt);

    // read boundary type
    typename Indices::BoundaryFormulation bcType;
    problem().boundaryFormulation(bcType, *isIt);
    if (bcType == Indices::saturation)
    {
        Scalar satBound = primaryVariablesOnBoundary[Indices::contiWEqIdx];
        if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
        {
            Scalar pcBound = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPosFace, *eIt),
                    satBound);
            switch (pressureType)
            {
            case pw:
            {
                pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx] + pcBound;
                break;
            }
            case pn:
            {
                pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx] - pcBound;
                pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                break;
            }
            }
        }
        else // capillarity neglected
            pressBound[wPhaseIdx] = pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];


        BCfluidState.satFlash(satBound, pressBound, problem().spatialParameters().porosity(globalPosFace, *eIt),
                            problem().temperatureAtPos(globalPosFace));

    }
    else if (bcType == Indices::concentration)
    {
        // saturation and hence pc and hence corresponding pressure unknown
        pressBound[wPhaseIdx] = pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
        Scalar Z1Bound = primaryVariablesOnBoundary[Indices::contiWEqIdx];
        BCfluidState.update(Z1Bound, pressBound, problem().spatialParameters().porosity(globalPosFace, *eIt),
                            problem().temperatureAtPos(globalPosFace));

        if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
        {
            Scalar pcBound = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPosFace, *eIt),
                    BCfluidState.saturation(wPhaseIdx));
            int maxiter = 3;
            //start iteration loop
            for(int iter=0; iter < maxiter; iter++)
            {
                //prepare pressures to enter flash calculation
                switch (pressureType)
                {
                case pw:
                {
                    pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                    pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx] + pcBound;
                    break;
                }
                case pn:
                {
                    pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx] - pcBound;
                    pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                    break;
                }
                }

                //store old pc
                Scalar oldPc = pcBound;
                //update with better pressures
                BCfluidState.update(Z1Bound, pressBound, problem().spatialParameters().porosity(globalPosFace, *eIt),
                                    problem().temperatureAtPos(globalPosFace));
                pcBound = MaterialLaw::pC(problem().spatialParameters().materialLawParams(globalPosFace, *eIt),
                        BCfluidState.saturation(wPhaseIdx));
                // TODO: get right criterion, do output for evaluation
                //converge criterion
                if (abs(oldPc-pcBound)<10)
                    iter = maxiter;
            }
        }
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Boundary Formulation neither Concentration nor Saturation??");
}

}
#endif
