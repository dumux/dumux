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
 * \brief Manages the domain for the lenhard problem. The domain consists of a
 *        grid and the states for the given entities.
 *
 * This has been modeled to match as closely as possible to the one
 * described at:
 *
 * Sheta, Hussam: "Simulation von Mehrphasenvorgaengen in poroesen
 *     Medien unter Einbeziehung von Hystereseeffekten", PhD theses,
 *     Braunschweig 1999, pp. 112
 */
#ifndef DUMUX_LENHARD_DOMAIN_HH
#define DUMUX_LENHARD_DOMAIN_HH

#include "lenhardstate.hh"

#include <dumux/auxiliary/basicdomain.hh>

#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/new_material/parkerlenhard.hh>

#include <dumux/material/properties.hh>

#include <dune/grid/sgrid.hh>
#include <dune/grid/onedgrid.hh>
//#include <dune/grid/yaspgrid.hh>
//#include <dumux/onedinndgrid/onedinndgrid.hh>
//#include <dumux/onedinndgrid/onedinndgrid.cc>

namespace Dune
{
namespace Lenhard
{
    typedef Dune::SGrid<1, 1> LenhardGrid;
//    typedef Dune::YaspGrid<1, 1> LenhardGrid;
//    typedef Dune::OneDGrid       LenhardGrid;

    /*!
     * \brief The domain for the Pw-Sn lenhard problem
     */
    template<class ScalarT>
    class PwSnLenhardDomain : public BasicDomain<LenhardGrid, ScalarT>
    {
        typedef PwSnLenhardDomain<ScalarT>                   ThisType;
        typedef BasicDomain<LenhardGrid, ScalarT>            ParentType;
        typedef ScalarT                                      Scalar;
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
            typedef LenhardGlobalState<ScalarT>                   GlobalState;
            typedef LenhardMediumState<VanGenuchtenState>         MediumState;
            typedef LenhardVertexState<MediumState>                 VertexState;
            typedef LenhardElementState<MediumState>                 ElementState;

            // The parker-lenhard hysteresis model
#if !defined USE_NODE_PARAMETERS
            typedef Dune::ParkerLenhard<ElementState, VanGenuchten>     ParkerLenhard;
#else // defined USE_NODE_PARAMETERS
            typedef Dune::ParkerLenhard<VertexState, VanGenuchten>   ParkerLenhard;
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

        typedef typename DomainTraits::LocalPosition    LocalPosition;
        typedef typename DomainTraits::GlobalPosition   GlobalPosition;
        typedef typename DomainTraits::FieldVector      FieldVector;
        typedef typename DomainTraits::FieldMatrix      FieldMatrix;

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
        PwSnLenhardDomain()
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                initGrid_();

                initGlobalState_();
                initMediaStates_();
                initElementStates_();
                initVertexStates_();
            };

        ~PwSnLenhardDomain()
            {
                delete coarseSand_;
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
#if defined USE_NODE_PARAMETERS
                return ParkerLenhard::pC(vertexState(globalVertIdx), Sw);
#else // !defined USE_NODE_PARAMETERS
#if defined USE_INTERFACE_CONDITION
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
/*                printf("krw(%f): %f, viscosity: %f\n",
                       Sw,
                       ParkerLenhard::krw(elementState(elementIdx), Sw),
                       viscosityW());
*/

#if defined USE_NODE_PARAMETERS
                return ParkerLenhard::krw(vertexState(globalVertIdx), Sw) / viscosityW();
#else // !defined USE_NODE_PARAMETERS
#ifdef USE_INTERFACE_CONDITION
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
/*                printf("krn(%f): %f, viscosity: %f\n",
                  1 - Sn,
                  ParkerLenhard::krn(elementState(elementIdx), 1 - Sn),
                  viscosityW());
*/

                // the capillary pressure models expect Sw as the
                // independend variable even for the non-wetting phase
                // mobility!
                Scalar Sw = 1 - Sn;

#if defined USE_NODE_PARAMETERS
                return ParkerLenhard::krn(vertexState(globalVertIdx), Sw) / viscosityN();
#else // !defined USE_NODE_PARAMETERS
#ifdef USE_INTERFACE_CONDITION
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
            { return vertStates_[ParentType::vertIdx(element, i)]; }
        const VertexState &vertexState(const Element &element, int i) const
            { return vertStates_[ParentType::vertIdx(element, i)]; }

        // given a vert, return it's state object
        VertexState &vertexState(const Vertex &vert)
            { return vertStates_[ParentType::vertIdx(vert)]; }
        const VertexState &vertexState(const Vertex &vert) const
            { return vertStates_[ParentType::vertIdx(vert)]; }

        const GlobalPosition &lowerLeft() const
            { return gridLowerLeft_; }
        const GlobalPosition &upperRight() const
            { return gridUpperRight_; }

        Scalar height() const
            { return gridUpperRight_[0] - gridLowerLeft_[0]; }

        bool onUpperBoundary(const GlobalPosition &pos) const
            { return pos[0] > gridUpperRight_[0] - eps_; }
        bool onLowerBoundary(const GlobalPosition &pos) const
            { return pos[0] < gridLowerLeft_[0] + eps_; }

    private:
        void initGrid_()
            {
                // specify the position and size of the grid in world
                // coordinates
                gridLowerLeft_[0] = 0;

                // column was 72 cm high, but we need some additional
                // space so that the boundary condition at the upper
                // boundary doesn't destroy the experiment
                gridUpperRight_[0] = 0.80;


                // the epsilon constant
                eps_ = 1e-8 * height();

                // create the grid
                Dune::FieldVector<int, dim> elementRes;
                elementRes[0] = 160;

                Grid *grid = new Grid(elementRes,
                                      gridLowerLeft_,
                                      gridUpperRight_);
                ParentType::setGrid(grid);
            }

        void initGlobalState_()
            {
                // initialize the globally uniform parameters
                globalState_ = new GlobalState();

                Water wettingFluid;
                Air   nonwettingFluid;
                globalState_->setDensityW(wettingFluid.density());
                globalState_->setDensityN(nonwettingFluid.density());

                globalState_->setViscosityW(wettingFluid.viscosity());
                globalState_->setViscosityN(nonwettingFluid.viscosity());

                FieldVector gravity(0);
                gravity[0] = -9.81;
                globalState_->setGravity(gravity);
            }

        void initMediaStates_()
            {
                // initializing the state of the porous medium
                coarseSand_ =  new MediumState();
                coarseSand_->setGlobalState(globalState_);
#if LENHARD_EXPERIMENT == 1
                coarseSand_->setSwr(0.0);
                coarseSand_->setSnr(0.25);
                coarseSand_->setPorosity(0.36);
                coarseSand_->setMicParams(VanGenuchtenState(0.00052, 4.63));
                coarseSand_->setMdcParams(VanGenuchtenState(0.00052, 4.63));
                coarseSand_->micParams().setVgMaxPC(coarseSand_->mdcParams().vgMaxPC());
                coarseSand_->setPermeability(3.37e-11);
#elif LENHARD_EXPERIMENT == 2
                coarseSand_->setSwr(0.17);
                coarseSand_->setSnr(0.25);
                coarseSand_->setPorosity(0.36);
                coarseSand_->setMdcParams(VanGenuchtenState(0.00042, 5.25));
//                coarseSand_->setMicParams(VanGenuchtenState(0.00042, 5.25));
                coarseSand_->setMicParams(VanGenuchtenState(0.00042*2, 5.25));
                coarseSand_->micParams().setVgMaxPC(coarseSand_->mdcParams().vgMaxPC());
                coarseSand_->setPermeability(3.37e-11);
#endif
            }

        void initElementStates_()
            {
                elementStates_.resize(ParentType::numElements());

                // initialize the element state objects depending on
                // wether the element is inside or outside of the lenhard of
                // fine sand.
                GlobalPosition elementCenter;
                ElementIterator it = ParentType::elementBegin();
                ElementIterator endit = ParentType::elementEnd();
                for (; it != endit; ++it) {
                    ParentType::elementCenter(*it, elementCenter);
                    int arrayPos = ParentType::elementIdx(*it);
                    elementStates_[arrayPos].setMediumState(coarseSand_);
                }

            }

        void initVertexStates_()
            {
                vertStates_.resize(ParentType::numVertices());

#if defined USE_NODE_PARAMETERS
                VertexIterator vertIt = ParentType::vertexBegin();
                const VertexIterator &endVertIt = ParentType::vertexEnd();
                for (; vertIt != endVertIt; ++vertIt) {
                    GlobalPosition pos;
                    ParentType::vertPosition(pos, *vertIt);
                    VertexState &vertState = vertexState(*vertIt);

                    vertState.setMediumState(coarseSand_);
                }

#else // ! defined USE_NODE_PARAMETERS
                // initialize the vert state objects
                /*
                // loop over all elements
                ElementIterator elementIt = ParentType::elementBegin();
                const ElementIterator &endElementIt = ParentType::elementEnd();
                for (; elementIt != endElementIt; ++elementIt) {
                    // loop over all faces of the current element
                    IntersectionIterator isIt = IntersectionIteratorGetter::begin(*elementIt);
                    const IntersectionIterator &endFaceIt = IntersectionIteratorGetter::end(*elementIt);
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
                */
#endif
            }

#if !defined USE_NODE_PARAMETERS
        //  mark all vertices on a face as as belonging to an
        //  inhomogenity
        /*
        void markInterfaceVertices_(IntersectionIterator &interisIt)
            {
                const Element &element = *interisIt.inside();
                const ReferenceElement &refElem =
                    DomainTraits::referenceElement(element.geometry().type());

                int inIdx = elementIdx(*interisIt.inside());
                int outIdx = elementIdx(*interisIt.outside());

                int faceIdx = interisIt.numberInSelf();
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
        */
#endif // !defined USE_NODE_PARAMETERS


        // the actual type which stores the element states.
        typedef std::vector<VertexState>  VertexStateArray_;

        // the actual type which stores the element states.
        typedef std::vector<ElementState>    ElementStateArray_;

        // the lower left and upper right coordinates of the complete
        // grid
        GlobalPosition gridLowerLeft_;
        GlobalPosition gridUpperRight_;

        // the lower left and upper right coordinates of the fine sand
        // lenhard
        GlobalPosition lenhardLowerLeft_;
        GlobalPosition lenhardUpperRight_;

        // global state and states of the media
        GlobalState *globalState_;
        MediumState *coarseSand_;

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
