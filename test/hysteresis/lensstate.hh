/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file LensState.hh The required state classes for the lens problem
 */
#ifndef DUMUX_LENSSTATE_HH
#define DUMUX_LENSSTATE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/new_material/statehelpermacros.hh>
#include <dumux/new_material/parkerlenhardstate.hh>
#include <dumux/new_material/parkerlenhard.hh>

#include <vector>

namespace Dune
{
namespace Lens
{
    // parameters which are globally constant
    template <class ScalarT>
    class LensGlobalState
    {
    public:
        typedef ScalarT                          Scalar;
        typedef Dune::FieldVector<Scalar, 2>     NVector;
        typedef Dune::FieldMatrix<Scalar, 2, 2>  NxNMatrix;

        // Density of the wetting phase
        PROPERTY(Scalar, densityW, setDensityW);

        // Density of the non-wetting phase
        PROPERTY(Scalar, densityN, setDensityN);

        // the gravity vector
        PROPERTY(NVector, gravity, setGravity);

        // the viscosity of the wetting fluid
        PROPERTY(Scalar, viscosityW, setViscosityW);

        // the viscosity of the non-wetting fluid
        PROPERTY(Scalar, viscosityN, setViscosityN);
    };

    // parameters which are constant for a porous medium
    template <class CapPressureParamsT>
    class LensMediumState
    {
    public:
        typedef CapPressureParamsT                  CapPressureParams;
        typedef typename CapPressureParams::Scalar  Scalar;
        typedef Dune::FieldVector<Scalar, 2>        NVector;
        typedef Dune::FieldMatrix<Scalar, 2, 2>     NxNMatrix;
        typedef LensGlobalState<Scalar>             GlobalState;


        LensMediumState()
            {
                Swr_ = Snr_ = 0;
                porosity_ = 0;
                globalState_ = NULL;
            }

        // Low level hysteresis model parameters (i.e. parameters
        // of the underlying model for the main imbibition and
        // main drainage curves)
        MUTABLE_PROPERTY(CapPressureParams, micParams, setMicParams);
        MUTABLE_PROPERTY(CapPressureParams, mdcParams, setMdcParams);

        // Twophase parameters
        PARAMETER(Scalar, Swr);
        void setSwr(Scalar Swr)
            {
                Swr_ = Swr;
                Snre_ = Snr_/(1 - Swr_);
            }
        PARAMETER(Scalar, Snr);
        void setSnr(Scalar Snr)
            {
                Snr_ = Snr;
                Snre_ = Snr_/(1 - Swr_);
            }
        PARAMETER(Scalar, Snre);

        // absolute permeability
        PROPERTY(NxNMatrix, permeability, setPermeability);
        void setPermeability(Scalar permeability)
            {
                permeability_[0][0] = permeability_[1][1] = permeability;
                permeability_[0][1] = permeability_[1][0] = 0;
            }

        // the absolute porosity of the medium
        PROPERTY(Scalar, porosity, setPorosity);

        // pointer to the global properties
        PTR_PROPERTY(GlobalState, globalState, setGlobalState);
    };


#if !USE_NODE_PARAMETERS
    // vert dependend parameters for element centered parameter storage
    template <class MediumStateT>
    class LensVertexState
    {
    public:
        LensVertexState()
            {
                neighbours_ = NULL;
            };

        ~LensVertexState()
            {
                delete neighbours_;
            };

        bool isOnInterface() const
            { return neighbours_ != NULL; }

        void addNeighbourElementIdx(int elementIdx)
            {
                if (neighbours_ == NULL) {
                    // allocate the array for the neighbouring element
                    // indices.  in order not having to reallocate the
                    // array every time a new index is added, we
                    // reserve 4 slots from the beginning.
                    neighbours_ = new NeighbourIdxArray_;
                    neighbours_->reserve(4);
                }

                // make sure the elementIdx is not already in the array
                for (unsigned i = 0; i < neighbours_->size(); ++i)
                    if ((*neighbours_)[i] == elementIdx)
                        return;

                // append the element index
                neighbours_->push_back(elementIdx);
            }

        int numNeighbourElements() const
            {
                return neighbours_ == NULL
                    ? 0
                    : neighbours_->size();
            }

        int neighbourElementIdx(int localIdx) const
            {
                assert(neighbours_ != NULL);
                return neighbours_->at(localIdx);
            }

    private:
        typedef std::vector<int> NeighbourIdxArray_;
        NeighbourIdxArray_ *neighbours_;
    };
#else // USE_NODE_PARAMETERS
    template <class MediumStateT>
    class LensElementState
    {
    public:
        typedef MediumStateT                              MediumState;

        typedef typename MediumState::Scalar              Scalar;
        typedef typename MediumState::NVector             NVector;
        typedef typename MediumState::NxNMatrix           NxNMatrix;

        // pointer to the medium uniform properties
        const MediumState *mediumState() const
            { return mediumState_; }
        void setMediumState(const MediumState *ms)
            { mediumState_ = ms; }

        //////////////////////////////////
        // medium uniform parameters (implemented as proxy parameters)
        PROXY_PARAMETER(mediumState(), NxNMatrix, permeability);
        PROXY_PARAMETER(mediumState(), Scalar, porosity);
        //////////////////////////////////

    private:
        const MediumState *mediumState_;
    };
#endif


    // element/vert dependent parameters
#if USE_NODE_PARAMETERS
#define MAIN_STATE_CLASS LensVertexState
#else
#define MAIN_STATE_CLASS LensElementState
#endif

    template <class MediumStateT>
    class MAIN_STATE_CLASS
    {
    public:
        typedef MediumStateT                              MediumState;
        typedef typename MediumState::GlobalState         GlobalState;

        typedef typename MediumState::Scalar              Scalar;
        typedef typename MediumState::NVector             NVector;
        typedef typename MediumState::NxNMatrix           NxNMatrix;
        typedef typename MediumState::CapPressureParams   CapPressureParams;
        typedef Dune::PLScanningCurve<Scalar>             ScanningCurve;

        MAIN_STATE_CLASS()
        {
            init_();
        };

        MAIN_STATE_CLASS(const MAIN_STATE_CLASS &s)
        {
            // we don't really copy the other element state, but only
            // use the same medium and global parameters.  (The sole
            // purpose of the copy constructor is to make STL
            // containers work properly.)
            init_(s.mediumState_);
        }

        MAIN_STATE_CLASS(const MediumState *mediumState)
        {
            init_(mediumState);
        }


        ~MAIN_STATE_CLASS()
            { }


        // pointer to the medium uniform properties
        const MediumState *mediumState() const
            { return mediumState_; }
        void setMediumState(const MediumState *ms)
            { mediumState_ = ms; }

        //////////////////////////////////
        // medium uniform parameters (implemented as proxy parameters)
        PROXY_PARAMETER(mediumState(), NxNMatrix, permeability);
        PROXY_PARAMETER(mediumState(), Scalar, porosity);
        //////////////////////////////////

        //////////////////////////////////
        // medium uniform parameters (implemented as proxy parameters)
        PROXY_PARAMETER(mediumState(), CapPressureParams, micParams);
        PROXY_PARAMETER(mediumState(), CapPressureParams, mdcParams);
        PROXY_PARAMETER(mediumState(), Scalar, Swr);
        PROXY_PARAMETER(mediumState(), Scalar, Snr);
        PROXY_PARAMETER(mediumState(), Scalar, Snre);
        //////////////////////////////////

        //////////////////////////////////
        // element specific properties
        MUTABLE_PTR_PROPERTY(ScanningCurve, mdc, setMdc); // main drainage scanning curve
        MUTABLE_PTR_PROPERTY(ScanningCurve, pisc, setPisc); // primary imbibition scanning curve
        MUTABLE_PTR_PROPERTY(ScanningCurve, csc, setCsc); // current scanning curve
        PROPERTY(Scalar, Snrei, setSnrei); // current effective residual non-wetting saturation
        //////////////////////////////////

    private:
        void init_(const MediumState *mediumState = NULL)
            {
                mediumState_ = mediumState;

                Snrei_ = 0;
                mdc_ = new ScanningCurve();
                pisc_ = csc_ = NULL;
            }

        const MediumState *mediumState_;
    };
}
}

#endif
