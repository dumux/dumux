/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
/*!
 * \file
 * \brief Provides the required interfaces for 2p2c problems.
 */
#ifndef DUMUX_BASIC_DOMAIN_HH
#define DUMUX_BASIC_DOMAIN_HH

#include <dumux/auxiliary/basicdomain.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>


namespace Dune
{
/*!
 * \brief Provides the required interfaces for a 2p2c problems.
 *
 * The Domains for the specific models are derived from this
 * class.
 */
template<class ModelTraitsT,
         class Implementation,
         class GridT,
         class ScalarT>
class TwoPTwoCBaseDomain : public BasicDomain<GridT, ScalarT>
{
    typedef BasicDomain<GridT, ScalarT> ParentT;

public:
    typedef typename ParentT::DomainTraits         DomainTraits;
    typedef typename ModelTraitsT                  ModelTraits;
    typedef Dune::P1BoxTraits<ScalarT,
                              GridT,
                              ModelTraitsT::numEq> BoxTraits;


private:
    // some types from the traits for convenience (The
    // DomainTraits are inherited from the BasicDomain)
    typedef typename DomainTraits::Element               Element;
    typedef typename DomainTraits::IntersectionIterator  IntersectionIterator;
    typedef typename DomainTraits::LocalPosition         LocalPosition;
    typedef typename DomainTraits::GlobalPosition        GlobalPosition;
    typedef typename BoxTraits::FVElementGeometry        FVElementGeometry;

public:
    TwoPTwoCBaseDomain()
        : ParentT()
    {
    };

    TwoPTwoCBaseDomain(Grid *grid)
        : ParentT(grid)
    {
    };

    /*!
     * \brief The capillary pressure at sub-control-volume scvIdx.
     */
    template <class Variables>
    Scalar pC(const Variables         &vars,
              const Element           &element,
              const FVElementGeometry &fvElemGeom,
              int                      scvIdx) const
    {
        return Implementation::materialLaw().pC(vars.saturation(wPhase),
                                                fvElemGeom.scv[scvIdx].global,
                                                element,
                                                fvElemGeom.scv[scvIdx].local);
    };

    /*!
     * \brief The mass fraction of non-wetting component in wetting
     *        phase ("air in water")
     */
    template <class Variables>
    Scalar xAW(const Variables         &vars,
               const Element           &element,
               const FVElementGeometry &fvElemGeom,
               int                      scvIdx) const
    {
        return Implementation::multicomp().xWN(vars.pressure(nPhase),
                                               vars.temperature(nPhase));
    };

    /*!
     * \brief The mass fraction of wetting component in non-wetting
     *        phase ("water in non-wetting")
     */
    template <class Variables>
    Scalar xWN(const Variables       &vars,
               const Element           &element,
               const FVElementGeometry &fvElemGeom,
               int                      scvIdx) const
    {
        return Implementation::multicomp().xWN(vars.pressure(nPhase),
                                               vars.temperature(nPhase));
    };

    /*!
     * \brief The density of a phase.
     */
    template <class Variables>
    Scalar density(int phaseIdx,
                   const Variables         &vars,
                   const Element           &element,
                   const FVElementGeometry &fvElemGeom,
                   int                      scvIdx) const
    {
        if (phaseIdx == wPhase) {
            return
                Implementation::wettingPhase().density(vars.temperature(nPhase),
                                                       vars.pressure(wPhase),
                                                       vertDat.massFrac(nComp, wPhase));
        }
        else if (phaseIdx == nPhase) {
            return
                Implementation::nonwettingPhase().density(vars.temperature(nPhase),
                                                          vars.pressure(nPhase),
                                                          vertDat.massFrac(wComp, nPhase));
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Phase index "
                       << phaseIdx
                       << " is invalid for 2p2c models."););
};

/*!
 * \brief The phase mobility.
 */
template <class Variables>
Scalar mobility(int phaseIdx,
                const Variables         &vars,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
{
    if (phaseIdx == wPhase) {
        return
            Implementation::materialLaw().mobW(vars.saturation(wPhase),
                                               fvElemGeom.scv[scvIdx].global,
                                               element,
                                               fvElemGeom.scv[scvIdx].local,
                                               vars.temperature(wPhase),
                                               vars.pressure(wPhase));
    }
    else if (phaseIdx == nPhase) {
        Implementation::materialLaw().mobW(vars.saturation(nPhase),
                                           fvElemGeom.scv[scvIdx].global,
                                           element,
                                           fvElemGeom.scv[scvIdx].local,
                                           vars.temperature(nPhase),
                                           vars.pressure(nPhase));
    }
    else
        DUNE_THROW(Dune::InvalidStateException,
                   "Phase index "
                   << phaseIdx
                   << " is invalid for 2p2c models."););

};
};
}

#endif
