/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
 * \file
 * \brief Manages the domain for the lens problem. The domain consists of a
 *        grid and the states for the given entities.
 */
#ifndef DUMUX_LENS_DOMAIN_HH
#define DUMUX_LENS_DOMAIN_HH

#include "lensstate.hh"

#include <dumux/auxiliary/basicdomain.hh>

#include <dumux/new_models/2p/pwsnboxmodel.hh>

#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/new_material/parkerlenhard.hh>

#include <dumux/material/properties.hh>

#include <dune/grid/sgrid.hh>

namespace Dune
{
namespace Lens
{
typedef Dune::SGrid<2, 2> LensGrid;
//      typedef Dune::YaspGrid<2, 2> LensGrid;
//      typedef Dune::UGGrid<2> LensGrid;

/*!
 * \brief The domain for the Pw-Sn lens problem
 */
template<class ScalarT>
class PwSnLensDomain : public BasicDomain<LensGrid, ScalarT>
{
    typedef PwSnLensDomain<ScalarT>                               ThisType;
    typedef BasicDomain<LensGrid, ScalarT>                        ParentType;
    typedef ScalarT                                               Scalar;

public:
    // traits for the domain
    typedef typename ParentType::DomainTraits DomainTraits;

    // traits of the material relations
    struct MaterialTraits
    {
        // The basic capillary pressure model used for the
        // hystersis model.
        typedef RegularizedVanGenuchtenState<ScalarT>      VanGenuchtenState;
        typedef RegularizedVanGenuchten<VanGenuchtenState> VanGenuchten;

        // States for globally constant parameters, parameters
        // specfic for a medium, and element specific parameters
        typedef LensGlobalState<ScalarT>                   GlobalState;
        typedef LensMediumState<VanGenuchtenState>         MediumState;
        typedef LensVertexState<MediumState>                 VertexState;
        typedef LensElementState<MediumState>                 ElementState;

        // The parker-lenhard hysteresis model
#if !USE_NODE_PARAMETERS
        typedef Dune::ParkerLenhard<ElementState, VanGenuchten, USE_SPLINES>     ParkerLenhard;
#else // USE_NODE_PARAMETERS
        typedef Dune::ParkerLenhard<VertexState, VanGenuchten, USE_SPLINES>   ParkerLenhard;
#endif
    };


private:

    typedef typename DomainTraits::Grid                  Grid;
    typedef typename DomainTraits::Element               Element;
    typedef typename DomainTraits::ReferenceElement      ReferenceElement;

    typedef typename DomainTraits::Vertex                Vertex;

    typedef typename DomainTraits::ElementIterator       ElementIterator;
    typedef typename DomainTraits::VertexIterator        VertexIterator;

    typedef typename DomainTraits::IntersectionIterator  IntersectionIterator;

    typedef typename DomainTraits::LocalPosition         LocalPosition;
    typedef typename DomainTraits::GlobalPosition        GlobalPosition;
    typedef typename DomainTraits::FieldVector           FieldVector;
    typedef typename DomainTraits::FieldMatrix           FieldMatrix;

    typedef typename MaterialTraits::GlobalState   GlobalState;
    typedef typename MaterialTraits::MediumState   MediumState;
    typedef typename MaterialTraits::ElementState  ElementState;
    typedef typename MaterialTraits::VertexState   VertexState;

    typedef typename MaterialTraits::VanGenuchtenState VanGenuchtenState;
    typedef typename MaterialTraits::VanGenuchten      VanGenuchten;
    typedef typename MaterialTraits::ParkerLenhard     ParkerLenhard;

    enum {
        dim = DomainTraits::dim
    };

public:
    PwSnLensDomain()
    {
        Api::require<Api::BasicDomainTraits, DomainTraits>();

        initGrid_();

        initGlobalState_();
        initMediaStates_();
        initElementStates_();
        initVertexStates_();
    };

    ~PwSnLensDomain()
    {
        delete outerMedium_;
        delete lensMedium_;
        delete globalState_;
    }

    //! returns the body force vector within a element
    const FieldVector &gravity(const Element &element,
                               int elementIdx) const
    {
        return globalState_->gravity();
    }

    const FieldVector &gravity() const
    {
        return globalState_->gravity();
    }

    //! return the permeability tensor for a element
    void applyPermeabilityTensor(FieldVector &dest, const Element &element, const FieldVector &pressureGradient) const
    {
        elementState(elementIdx(element)).permeability().umv(pressureGradient, dest);
    }

    //! Return the capillary pressure for a given vert of a element
    Scalar pC(const Element &element,
              int elementIdx,
              int localVertIdx,
              int globalVertIdx,
              Scalar Sw) const
    {
#if USE_NODE_PARAMETERS
        return ParkerLenhard::pC(vertexState(globalVertIdx), Sw);
#else // !USE_NODE_PARAMETERS
#if USE_INTERFACE_CONDITION
        const VertexState &vertState = vertexState(globalVertIdx);
        if (!vertState.isOnInterface())
            return ParkerLenhard::pC(elementState(elementIdx),
                                     Sw);

        // find the minimum capilarry pressure around the
        // vert if the vert is on an interface
        Scalar minPc = 1e100;
        int nNeighbors = vertState.numNeighbourElements();
        for (int i = 0; i < nNeighbors; ++i) {
            const ElementState &ncState = elementState(vertState.neighbourElementIdx(i));
            minPc = std::min(minPc,
                             ParkerLenhard::pC(ncState, Sw));
        }

        return minPc;
#else // ! USE_INTERFACE_CONDITION
        return ParkerLenhard::pC(elementState(elementIdx),
                                 Sw);
#endif
#endif
    }

    // return the derivative of the capillary pressure regarding
    // the wetting phase saturation
    Scalar dpC_dSw(const Element &element,
                   int elementIdx,
                   int localVertIdx,
                   int globalVertIdx,
                   Scalar Sw) const
    {
        /* TODO
           const ElementState &state = elementState(elementIdx);
           Scalar a = ParkerLenhard::pC(state, Sw + 5e-4);
           Scalar b = ParkerLenhard::pC(state, Sw - 5e-4);
           return (a - b) / 1e-3;
        */
        return 0;
    }

    // return the capillary pressure for a given element
    Scalar porosity(const Element &element) const
    {
        return elementState(element).porosity();
    }

    // return the density of the wetting phase
    Scalar densityW() const
    {
        return globalState_->densityW();
    }

    // return the density of the non-wetting phase
    Scalar densityN() const
    {
        return globalState_->densityN();
    }

    // return the density of a phase given by an index. (lower
    // indices means that the phase is more wetting)
    Scalar density(int phase) const
    {
        return (phase == 0)? densityW() : densityN();
    }

    // return the viscosity of the wetting phase
    Scalar viscosityW() const
    {
        return globalState_->viscosityW();
    }

    // return the viscosity of the non-wetting phase
    Scalar viscosityN() const
    {
        return globalState_->viscosityN();
    }

    // return the viscosity of a phase given by an index. (lower
    // indices means that the phase is more wetting)
    Scalar viscosity(int phase) const
    {
        return (phase == 0)? viscosityW() : viscosityN();
    }

    // return the mobility of the wetting phase at a vert
    Scalar mobilityW(const Element &element,
                     int elementIdx,
                     int localVertIdx,
                     int globalVertIdx,
                     Scalar Sw) const
    {
#if USE_NODE_PARAMETERS
        return ParkerLenhard::krw(vertexState(globalVertIdx), Sw) / viscosityW();
#else // !USE_NODE_PARAMETERS
#if USE_INTERFACE_CONDITION
        const VertexState &vertState = vertexState(globalVertIdx);
        if (!vertState.isOnInterface())
            return ParkerLenhard::krw(elementState(elementIdx), Sw) / viscosityW();

        // find the minimum wetting phase mobility around the
        // vert if the vert is on an interface
        Scalar minMob = 1e100;
        int nNeighbors = vertState.numNeighbourElements();
        for (int i = 0; i < nNeighbors; ++i) {
            const ElementState &ncState = elementState(vertState.neighbourElementIdx(i));
            minMob = std::min(minMob,
                              ParkerLenhard::krw(ncState, Sw) / viscosityW());
        };

        return minMob;
#else // ! USE_INTERFACE_CONDITION
        return ParkerLenhard::krw(elementState(elementIdx), Sw) / viscosityW();
#endif
#endif
    }

    // return the mobility of the non-wetting phase at a vert
    Scalar mobilityN(const Element &element,
                     int elementIdx,
                     int localVertIdx,
                     int globalVertIdx,
                     Scalar Sn) const
    {
        // the capillary pressure models expect Sw as the
        // independend variable even for the non-wetting phase
        // mobility!
        Scalar Sw = 1 - Sn;

#if USE_NODE_PARAMETERS
        return ParkerLenhard::krn(vertexState(globalVertIdx), Sw) / viscosityN();
#else // !USE_NODE_PARAMETERS
#if USE_INTERFACE_CONDITION
        const VertexState &vertState = vertexState(globalVertIdx);
        if (!vertState.isOnInterface())
            return ParkerLenhard::krn(elementState(elementIdx), Sw) / viscosityN();

        // find the minimum wetting phase mobility around the
        // vert if the vert is on an interface
        Scalar minMob = 1e100;
        int nNeighbors = vertState.numNeighbourElements();
        for (int i = 0; i < nNeighbors; ++i) {
            const ElementState &ncState = elementState(vertState.neighbourElementIdx(i));
            minMob = std::min(minMob,
                              ParkerLenhard::krn(ncState, Sw) / viscosityN());
        };

        return minMob;
#else // ! USE_INTERFACE_CONDITION
        return ParkerLenhard::krn(elementState(elementIdx), Sw) / viscosityN();
#endif
#endif
    }

    // given a element index, return the corresponding element state
    ElementState &elementState(int index)
    { return elementStates_[index]; }
    const ElementState &elementState(int index) const
    { return elementStates_[index]; }

    ElementState &elementState(const Element &element)
    { return elementStates_[elementIdx(element)]; }
    const ElementState &elementState(const Element &element) const
    { return elementStates_[elementIdx(element)]; }

    // given a global vert index, return the corresponding state
    VertexState &vertexState(int index)
    { return vertStates_[index]; }
    const VertexState &vertexState(int index) const
    { return vertStates_[index]; }

    // given a element and a local vert index, return the corresponding state
    VertexState &vertexState(const Element &element, int i)
    { return vertStates_[ParentType::vertexIdx(element, i)]; }
    const VertexState &vertexState(const Element &element, int i) const
    { return vertStates_[ParentType::vertexIdx(element, i)]; }

    // given a vert, return it's state object
    VertexState &vertexState(const Vertex &vert)
    { return vertStates_[ParentType::vertexIdx(vert)]; }
    const VertexState &vertexState(const Vertex &vert) const
    { return vertStates_[ParentType::vertexIdx(vert)]; }

    const GlobalPosition &lowerLeft() const
    { return gridLowerLeft_; }
    const GlobalPosition &upperRight() const
    { return gridUpperRight_; }

    Scalar width() const
    { return gridUpperRight_[0] - gridLowerLeft_[0]; }
    Scalar height() const
    { return gridUpperRight_[1] - gridLowerLeft_[1]; }

    bool onUpperBoundary(const GlobalPosition &pos) const
    { return pos[1] > gridUpperRight_[1] - eps_; }
    bool onLowerBoundary(const GlobalPosition &pos) const
    { return pos[1] < gridLowerLeft_[1] + eps_; }
    bool onLeftBoundary(const GlobalPosition &pos) const
    { return pos[0] < gridLowerLeft_[0] + eps_; }
    bool onRightBoundary(const GlobalPosition &pos) const
    { return pos[0] > gridUpperRight_[0] - eps_; }

    // returns true iff a world coordinate is within the fine
    // sand lens.
    bool isInLens(const GlobalPosition &coord)
    {
        for (int i = 0; i < DomainTraits::dimWorld; ++i) {
            if (lensLowerLeft_[i] > coord[i] ||
                lensUpperRight_[i] < coord[i])
            {
                return false;
            }
        }
        return true;
    }

    /*
    // returns true iff the element corrosponding to the state is
    // within the fine sand lens.
    bool isInLens(const VertexState &state)
    {
    return state.mediumState() == lensMedium_;
    }
    */

private:
    void initGrid_()
    {
        // specify the position and size of the grid in world
        // coordinates
        gridLowerLeft_[0] = 0;
        gridLowerLeft_[1] = 0;
        gridUpperRight_[0] = 6;
        gridUpperRight_[1] = 4;

        // the epsilon constant
        eps_ = 1e-8 * width();

        // create the grid
        Dune::FieldVector<int, dim> elementRes;
#if USE_ORIG_PROB
        elementRes[0] = 48;
        elementRes[1] = 32;
#else
        elementRes[0] = CELLRES_X;
        elementRes[1] = CELLRES_Y;
#endif

        Grid *grid = new Grid(elementRes,
                              gridLowerLeft_,
                              gridUpperRight_);
        ParentType::setGrid(grid);

        // specify the location and the size of the fine sand
        // lens in world coordinates
        lensLowerLeft_[0] = width()*(1./6);
        lensLowerLeft_[1] = height()*(2./4);
        lensUpperRight_[0] = width()*(4./6);
        lensUpperRight_[1] = height()*(3./4);
        // make sure the global coordniates of the lens are
        // aligned to the grid's element boundaries in order to
        // make sure the interface occures at the same
        // position for sized lenses for vert and element based
        // material parameter storage. (this is not strictly
        // necessary, but makes it easier to compare the
        // different approaches.)
        lensLowerLeft_[0] = width()/elementRes[0]
            * (int) (lensLowerLeft_[0]*elementRes[0]/width())
            + (0.5 - eps_)/elementRes[0];
        lensLowerLeft_[1] = height()/elementRes[1]
            * (int) (lensLowerLeft_[1]*elementRes[1]/height())
            + (0.5 - eps_)/elementRes[1];
        lensUpperRight_[0] = width()/elementRes[0]
            * (int) (lensUpperRight_[0]*elementRes[0]/width())
            + (-0.5 + eps_)/elementRes[0];
        lensUpperRight_[1] = height()/elementRes[1]
            * (int) (lensUpperRight_[1]*elementRes[1]/height())
            + (-0.5 + eps_)/elementRes[1];
    }

    void initGlobalState_()
    {
        // initialize the globally uniform parameters
        globalState_ = new GlobalState();

        Water water;
        DNAPL dnapl;
        globalState_->setDensityW(water.density());
        globalState_->setDensityN(dnapl.density());

        globalState_->setViscosityW(water.viscosity());
        globalState_->setViscosityN(dnapl.viscosity());
        FieldVector gravity(0); gravity[1] = -9.81;
        globalState_->setGravity(gravity);
    }

    void initMediaStates_()
    {
        // initializing the state of the lens' medium
        lensMedium_ =  new MediumState();
        lensMedium_->setGlobalState(globalState_);
        lensMedium_->setSwr(0.18);
#if USE_ORIG_PROB
        lensMedium_->setSnr(0.0);
#else
        lensMedium_->setSnr(0.25);
#endif
        lensMedium_->setPorosity(0.4);
        lensMedium_->setMdcParams(VanGenuchtenState(0.00045, 7.3));
#if USE_DIFFERENT_MAIN_CURVES
        lensMedium_->setMicParams(VanGenuchtenState(0.00045*2, 7.5));
#else
        lensMedium_->setMicParams(VanGenuchtenState(0.00045, 7.3));
#endif
        lensMedium_->micParams().setVgMaxPC(lensMedium_->mdcParams().vgMaxPC());
        lensMedium_->setPermeability(9.05e-13);

        // initializing the state of the outer medium
        outerMedium_ =  new MediumState();
        outerMedium_->setGlobalState(globalState_);
        outerMedium_->setSwr(0.05);
#if USE_ORIG_PROB
        outerMedium_->setSnr(0.0);
#else
        outerMedium_->setSnr(0.20);
#endif
        outerMedium_->setPorosity(0.4);
        outerMedium_->setMdcParams(VanGenuchtenState(0.0037, 4.7));
#if USE_DIFFERENT_MAIN_CURVES
        outerMedium_->setMicParams(VanGenuchtenState(0.0037*2, 4.8));
#else
        outerMedium_->setMicParams(VanGenuchtenState(0.0037, 4.7));
#endif

        outerMedium_->micParams().setVgMaxPC(outerMedium_->mdcParams().vgMaxPC());
        outerMedium_->setPermeability(4.6e-10);

#if 0
        lensMedium_->setSwr(outerMedium_->Swr());
        lensMedium_->setSnr(outerMedium_->Snr());
        lensMedium_->setPorosity(outerMedium_->porosity());
        lensMedium_->setMicParams(outerMedium_->micParams());
        lensMedium_->setMdcParams(outerMedium_->mdcParams());
        lensMedium_->setPermeability(outerMedium_->permeability());
#endif

    }

    void initElementStates_()
    {
        elementStates_.resize(ParentType::numElements());

        // initialize the element state objects depending on
        // wether the element is inside or outside of the lens of
        // fine sand.
        GlobalPosition elementCenter;
        ElementIterator it = ParentType::elementBegin();
        ElementIterator endit = ParentType::elementEnd();
        for (; it != endit; ++it) {
            ParentType::elementCenter(*it, elementCenter);
            int arrayPos = ParentType::elementIdx(*it);
            if (isInLens(elementCenter))
                elementStates_[arrayPos].setMediumState(lensMedium_);
            else
                elementStates_[arrayPos].setMediumState(outerMedium_);
        }
    }

    void initVertexStates_()
    {
        vertStates_.resize(ParentType::numVertices());
#if USE_NODE_PARAMETERS
        VertexIterator vertIt = ParentType::vertexBegin();
        const VertexIterator &endVertIt = ParentType::vertexEnd();
        for (; vertIt != endVertIt; ++vertIt) {
            GlobalPosition pos;
            ParentType::vertexPosition(pos, *vertIt);
            VertexState &vertState = vertexState(*vertIt);
            if (isInLens(pos))
                vertState.setMediumState(lensMedium_);
            else
                vertState.setMediumState(outerMedium_);
        }

#else // !USE_NODE_PARAMETERS
        // initialize the vert state objects

        // loop over all elements
        ElementIterator elementIt = ParentType::elementBegin();
        const ElementIterator &endElementIt = ParentType::elementEnd();
        for (; elementIt != endElementIt; ++elementIt) {
            // loop over all faces of the current element
            IntersectionIterator isIt = elementIt->ileafbegin();
            const IntersectionIterator &endFaceIt = elementIt->ileafend();
            for (; isIt != endFaceIt; ++ isIt) {
                // make sure the face not on the boundary of
                // the grid
                if (!isIt.neighbor())
                    continue;

                // check whether both elements share the same
                // medium. if yes we don't need to do anything
                // about the current face.
                const ElementState &inState = elementState(*isIt.inside());
                const ElementState &outState = elementState(*isIt.outside());
                if (inState.mediumState() != outState.mediumState())
                    continue;

                // alright, the face is on a inhomogenity, so
                // we have to mark all its vertices
                markInterfaceVertices_(isIt);
            }
        }
#endif
    }

#if !USE_NODE_PARAMETERS
    //  mark all vertices on a face as as belonging to an
    //  inhomogenity
    void markInterfaceVertices_(IntersectionIterator &interisIt)
    {
        const Element &element = *interisIt.inside();
        const ReferenceElement &refElem =
            DomainTraits::referenceElement(element.geometry().type());

        int inIdx = elementIdx(*interisIt.inside());
        int outIdx = elementIdx(*interisIt.outside());

        int faceIdx = interisIt->indexInInside();
        int numVerticesOfFace = refElem.size(faceIdx, 1, dim);
        for (int vertInFace = 0;
             vertInFace < numVerticesOfFace;
             vertInFace++)
        {
            int vertIdxInElement = refElem.subEntity(faceIdx, 1, vertInFace, dim);

            vertexState(*interisIt.inside(), vertIdxInElement).addNeighbourElementIdx(inIdx);
            vertexState(*interisIt.inside(), vertIdxInElement).addNeighbourElementIdx(outIdx);
        }
    }
#endif // !USE_NODE_PARAMETERS


    // the actual type which stores the element states.
    typedef std::vector<VertexState>  VertexStateArray_;

    // the actual type which stores the element states.
    typedef std::vector<ElementState>    ElementStateArray_;

    // the lower left and upper right coordinates of the complete
    // grid
    GlobalPosition gridLowerLeft_;
    GlobalPosition gridUpperRight_;

    // the lower left and upper right coordinates of the fine sand
    // lens
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    // global state and states of the media
    GlobalState *globalState_;
    MediumState *outerMedium_;
    MediumState *lensMedium_;

    // stores a state for each vert
    VertexStateArray_ vertStates_;

    // stores a state for each vert
    ElementStateArray_  elementStates_;

    // A small epsilon value and the current simulated
    // time. FIXME/TODO: should probably not be here..
    Scalar eps_;
};
}
}


#endif
