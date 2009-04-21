/*****************************************************************************
 *   Copyright (C) 2008,2009 by Melanie Darcis                               *
 *                              Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: melanie.darcis _at_ iws.uni-stuttgart.de                         *
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
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
#ifndef DUMUX_2P2CNI_VERTEX_DATA_HH
#define DUMUX_2P2CNI_VERTEX_DATA_HH

#include <dumux/new_models/2p2c/2p2cvertexdata.hh>

namespace Dune
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
template <class TwoPTwoCNITraits, 
          class Problem>
class TwoPTwoCNIVertexData : public TwoPTwoCVertexData<TwoPTwoCNITraits, Problem>
{
    typedef TwoPTwoCNITraits Tr;
    typedef typename Problem::DomainTraits::Scalar Scalar;
    typedef typename Problem::DomainTraits::Grid Grid;

    typedef typename Grid::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, Tr::numEq>      SolutionVector;
    typedef Dune::FieldVector<Scalar, Tr::numPhases>  PhasesVector;

    typedef Dune::FieldVector<Scalar, Grid::dimensionworld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, Grid::dimension>       LocalPosition;

    typedef TwoPTwoCVertexData<TwoPTwoCNITraits, Problem> ParentType;

    enum {
        dim = Grid::dimension
    };

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    template <class JacobianImp>
    void update(const SolutionVector   &sol,
                const Element          &element,
                int                     vertIdx,
                bool                    isOldSol,
                JacobianImp            &jac) 
    {
        // vertex update data for the mass balance
        ParentType::update(sol,
                           element,
                           vertIdx,
                           isOldSol,
                           jac);
        
        // data for the energy equation
        const LocalPosition &local =
            Problem::DomainTraits::referenceElement(element.type()).position(vertIdx,
                                                                             dim);
        const GlobalPosition &global =
            element.geometry().corner(vertIdx);

        temperature = sol[Tr::temperatureIdx];
        
        heatCond = jac.problem().soil().heatCond(global, 
                                                 element,
                                                 local,
                                                 this->saturation[Tr::wPhase]);
        
        enthalpy[Tr::wPhase] = jac.problem().wettingPhase().enthalpy(temperature,
                                                                     this->pressure[Tr::wPhase],
                                                                     this->massfrac[Tr::nComp][Tr::wPhase]);
        enthalpy[Tr::nPhase] = jac.problem().nonwettingPhase().enthalpy(temperature,
                                                                        this->pressure[Tr::nPhase],
                                                                        this->massfrac[Tr::wComp][Tr::nPhase]);
        intEnergy[Tr::wPhase] = jac.problem().wettingPhase().intEnergy(temperature,
                                                                       this->pressure[Tr::wPhase],
                                                                       this->massfrac[Tr::nComp][Tr::wPhase]);
        intEnergy[Tr::nPhase] = jac.problem().nonwettingPhase().intEnergy(temperature,
                                                                          this->pressure[Tr::nPhase],
                                                                          this->massfrac[Tr::wComp][Tr::nPhase]);
    }

    PhasesVector intEnergy; //!< Internal energy.
    PhasesVector enthalpy;  //!< Enthalpy.
    Scalar       temperature; //!< The temperature. We assume thermal equilibrium
    Scalar       heatCond; //!< Total heat conductivity.
};

} // end namepace

#endif
