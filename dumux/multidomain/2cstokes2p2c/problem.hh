// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
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
 * \brief The problem class for the coupling of a isothermal two-component Stokes
 *        and a isothermal two-phase two-component Darcy model.
 */

#ifndef DUMUX_2C_STOKES_2P2C_PROBLEM_HH
#define DUMUX_2C_STOKES_2P2C_PROBLEM_HH

#include <dumux/freeflow/boundarylayermodel.hh>
#include <dumux/freeflow/masstransfermodel.hh>
#include <dumux/freeflow/stokesnc/properties.hh>
#include <dumux/multidomain/problem.hh>
#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCStokesTwoCModel
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \brief The problem class for the coupling of a isothermal two-component Stokes
 *        and a isothermal two-phase two-component Darcy model.
 */
template <class TypeTag>
class TwoCStokesTwoPTwoCProblem : public MultiDomainProblem<TypeTag>
{
    typedef MultiDomainProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) Stokes2cTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TwoPTwoCTypeTag;
    typedef typename GET_PROP_TYPE(Stokes2cTypeTag, Indices) Stokes2cIndices;
    typedef typename GET_PROP_TYPE(TwoPTwoCTypeTag, Indices) TwoPTwoCIndices;

    // Stokes
    enum { numComponents1 = Stokes2cIndices::numComponents };
    enum { // component indices
        transportCompIdx1 = Stokes2cIndices::transportCompIdx, //!< Index of transported component
        phaseCompIdx1 = Stokes2cIndices::phaseCompIdx    //!< Index of main component of the phase
    };
    // Darcy
    enum { // phase indices
        wPhaseIdx2 = TwoPTwoCIndices::wPhaseIdx,         //!< Index for the liquid phase
    };

public:
    //! The constructor
    template<class GridView>
    TwoCStokesTwoPTwoCProblem(TimeManager &timeManager,
                              GridView gridView)
    : ParentType(timeManager, gridView)
    {
        blModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, BoundaryLayer, Model);
        massTransferModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, MassTransfer, Model);
    }

    /*!
     * \brief Returns a BoundaryLayerModel object
     *
     * \param cParams a parameter container
     * \param scvIdx1 The local index of the sub-control volume of the Stokes domain
     */
    template<typename CParams>
    BoundaryLayerModel<TypeTag> evalBoundaryLayerModel(CParams cParams, const int scvIdx1) const
    {
        // current position + additional virtual runup distance
        const Scalar distance = cParams.fvGeometry1.subContVol[scvIdx1].global[0]
                                + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, Offset);
        BoundaryLayerModel<TypeTag> boundaryLayerModel(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity),
                                                       distance,
                                                       cParams.elemVolVarsCur1[scvIdx1].kinematicViscosity(),
                                                       blModel_);
        if (blModel_ == 1)
            boundaryLayerModel.setConstThickness(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, ConstThickness));
        if (blModel_ >= 4)
            boundaryLayerModel.setYPlus(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, YPlus));
        if (blModel_ >= 5)
            boundaryLayerModel.setRoughnessLength(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, RoughnessLength));
        if (blModel_ == 7)
            boundaryLayerModel.setHydraulicDiameter(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryLayer, HydraulicDiameter));

        return boundaryLayerModel;
    }

    /*!
     * \brief Returns the concentration gradient through the boundary layer
     *
     * \param cParams a parameter container
     * \param scvIdx1 The local index of the sub-control volume of the Stokes domain
     */
    template<typename CParams>
    Scalar evalBoundaryLayerConcentrationGradient(CParams cParams, const int scvIdx1) const
    {
        static_assert(numComponents1 == 2,
                      "This coupling condition is only implemented for two components.");
        Scalar massFractionOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
        Scalar M1 = FluidSystem::molarMass(transportCompIdx1);
        Scalar M2 = FluidSystem::molarMass(phaseCompIdx1);
        Scalar X2 = 1.0 - massFractionOut;
        Scalar massToMoleDenominator = M2 + X2*(M1 - M2);
        Scalar moleFractionOut = massFractionOut * M2 /massToMoleDenominator;

        Scalar normalMoleFracGrad = cParams.elemVolVarsCur1[scvIdx1].moleFraction(transportCompIdx1)
                                    - moleFractionOut;
        return normalMoleFracGrad / asImp_().evalBoundaryLayerModel(cParams, scvIdx1).massBoundaryLayerThickness();
    }

    /*!
     * \brief Returns the mass transfer coefficient
     *
     * \param cParams a parameter container
     * \param scvIdx1 The local index of the sub-control volume of the Stokes domain
     * \param scvIdx2 The local index of the sub-control volume of the Darcy domain
     */
    template<typename CParams>
    Scalar evalMassTransferCoefficient(CParams cParams, const int scvIdx1, const int scvIdx2) const
    {
        MassTransferModel<TypeTag> massTransferModel(cParams.elemVolVarsCur2[scvIdx2].saturation(wPhaseIdx2),
                                                     cParams.elemVolVarsCur2[scvIdx2].porosity(),
                                                     asImp_().evalBoundaryLayerModel(cParams, scvIdx1).massBoundaryLayerThickness(),
                                                     massTransferModel_);
        if (massTransferModel_ == 1)
            massTransferModel.setMassTransferCoeff(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, MassTransfer, Coefficient));
        if (massTransferModel_ == 2 || massTransferModel_ == 4)
            massTransferModel.setCharPoreRadius(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, MassTransfer, CharPoreRadius));
        if (massTransferModel_ == 3)
            massTransferModel.setCapillaryPressure(cParams.elemVolVarsCur2[scvIdx2].capillaryPressure());

        return massTransferModel.massTransferCoefficient();
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    unsigned int blModel_;
    unsigned int massTransferModel_;
};

} // namespace Dumux

#endif // DUMUX_2C_STOKES_2P2C_PROBLEM_HH
