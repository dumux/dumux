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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCDarcysLaw
 */
#ifndef DUMUX_2P1C_SPURIOUS_FLUX_BLOCKING_DARCYS_LAW_HH
#define DUMUX_2P1C_SPURIOUS_FLUX_BLOCKING_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/darcyslaw.hh>

namespace Dumux
{

namespace Properties
{
     NEW_PROP_TAG(UseBlockingOfSpuriousFlow); //!< Determines whether blocking of spurious flow is used or not.
}

/*!
 * \ingroup TwoPOneCModel
 * \brief Specialization of Darcy's Law for the two-phase one-component model, including a the possibility <BR>
 *        to block spurious fluxes of cold water into the steam zone, which can improve the model's convergence
 *        behavior (Gudbjerg et al., 2005) \cite gudbjerg2004.
 */
template <class TypeTag>
class TwoPOneCDarcysLaw : public DarcysLaw<TypeTag>
{
    using ParentType = DarcysLaw<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElemFluxVarCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVarCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using CoordScalar = typename GridView::ctype;

    enum { dim = GridView::dimension};
    enum { dimWorld = GridView::dimensionworld};

    // copy some indices for convenience
    enum {
        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
    };

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    //! Compute the Darcy flux and use a blocking factor, if specified.
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const IndexType phaseIdx,
                       const ElemFluxVarCache& elemFluxVarCache)
    {
        const Scalar flux = ParentType::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarCache);

        // only block wetting-phase (i.e. liquid water) fluxes
        if((!GET_PROP_VALUE(TypeTag, UseBlockingOfSpuriousFlow)) || phaseIdx != wPhaseIdx)
            return flux;

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const Scalar factor = std::signbit(flux) ? factor_(outsideVolVars, insideVolVars, phaseIdx) :
                                                   factor_(insideVolVars, outsideVolVars, phaseIdx);

        return flux * factor;
    }

private:
    /*!
     * \brief Calculate the blocking factor which prevents spurious cold water fluxes into the steam zone (Gudbjerg et al., 2005) \cite gudbjerg2004 <BR>
     *
     * \param up The upstream volume variables
     * \param dn The downstream volume variables
     * \param phaseIdx The index of the fluid phase
     */
    static Scalar factor_(const VolumeVariables &up, const VolumeVariables &dn,const  int phaseIdx)
    {
        Scalar factor = 1.0;

        const Scalar tDn = dn.temperature(); //temperature of the downstream SCV (where the cold water is potentially intruding into a steam zone)
        const Scalar tUp = up.temperature(); //temperature of the upstream SCV

        const Scalar sgDn = dn.saturation(nPhaseIdx); //gas phase saturation of the downstream SCV
        const Scalar sgUp = up.saturation(nPhaseIdx); //gas phase saturation of the upstream SCV

        bool upIsNotSteam = false;
        bool downIsSteam = false;
        bool spuriousFlow = false;

        if(sgUp <= 1e-5)
            upIsNotSteam = true;

        if(sgDn > 1e-5)
            downIsSteam = true;

        if(upIsNotSteam && downIsSteam  && tDn > tUp && phaseIdx == wPhaseIdx)
          spuriousFlow = true;

        if(spuriousFlow)
        {
          Scalar deltaT = tDn - tUp;

          if((deltaT) > 15 )
            factor = 0.0 ;

          else
            factor = 1-(deltaT/15);
        }
        return factor;
    }
};

} // end namespace Dumux

#endif // DUMUX_2P1C_SPURIOUS_FLUX_BLOCKING_DARCYS_LAW_HH
