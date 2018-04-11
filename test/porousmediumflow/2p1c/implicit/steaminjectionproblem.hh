// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later vesion.                                      *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Non-isothermal steam injection test problem for the 2p1cni model.
 *
 */
#ifndef DUMUX_STEAM_INJECTIONPROBLEM_HH
#define DUMUX_STEAM_INJECTIONPROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/2p1c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/2p1c.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>

#include "steaminjectionspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblemTypeTag, INHERITS_FROM(TwoPOneCNI, InjectionProblemSpatialParams));
NEW_TYPE_TAG(TwoPOneCNIBoxTypeTag, INHERITS_FROM(BoxModel, InjectionProblemTypeTag));
NEW_TYPE_TAG(TwoPOneCNICCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, InjectionProblemTypeTag));

SET_TYPE_PROP(InjectionProblemTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InjectionProblemTypeTag, Problem, InjectionProblem<TypeTag>);


// Set fluid configuration
SET_PROP(InjectionProblemTypeTag, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using H2OType = Dumux::Components::TabulatedComponent<Dumux::Components::H2O<Scalar> >;
public:
    using type = Dumux::FluidSystems::TwoPOneC<Scalar, H2OType >;
};


//Define whether spurious cold-water flow into the steam is blocked
SET_BOOL_PROP(InjectionProblemTypeTag, UseBlockingOfSpuriousFlow, true);
}

/*!
 * \ingroup TwoPOneCTests
 * \brief Non-isothermal 2D problem where steam is injected on the lower left side of the domain.
 *
 * This problem uses the \ref TwoPOneC model.
 *
 *  */
template <class TypeTag>
class InjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    // copy some indices for convenience
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx,

        // phase state
        wPhaseOnly = Indices::wPhaseOnly
    };

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    /*!
     * \brief The constructor
     *
     */
    InjectionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        FluidSystem::init();
    }

    /*!
     * \name Problem parameters
     */
    // \{


    //! \copydoc Dumux::FVProblem::source()
    NumEqVector source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
          return NumEqVector(0.0);
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        if(globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_ || globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
           bcTypes.setAllDirichlet();
        else
           bcTypes.setAllNeumann();

         return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& ipGlobal = scvf.ipGlobal();

        if (ipGlobal[0] < eps_)
        {
            if(ipGlobal[1] > 2.0 - eps_ && ipGlobal[1] < 3.0 + eps_)
            {
                const Scalar massRate = 1e-1;
                values[conti0EqIdx] = -massRate;
                values[energyEqIdx] = -massRate * 2690e3;
            }
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        const Scalar densityW = 1000.0;
        values[pressureIdx] = 101300.0 + (this->fvGridGeometry().bBoxMax()[1] - globalPos[1])*densityW*9.81; // hydrostatic pressure
        values[switchIdx] = 283.13;

        values.setState(wPhaseOnly);

        return values;
    }

private:

    static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
