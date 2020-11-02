// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup TwoPModel
 * \copydoc Dumux::TwoPScvSaturationReconstruction
 */

#ifndef DUMUX_2P_SCV_SATURATION_RECONSTRUCTION_HH
#define DUMUX_2P_SCV_SATURATION_RECONSTRUCTION_HH

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Class that computes the nonwetting saturation in an scv from the saturation
 *        at the global degree of freedom.
 *
 * This is only necessary in conjunction with the box scheme where the degrees of
 * freedom lie on material interfaces. There the nonwetting phase saturation is
 * generally discontinuous.
 */
template<DiscretizationMethod M, bool enableReconstruction>
class TwoPScvSaturationReconstruction
{
public:
    /*!
     * \brief Compute the nonwetting phase saturation in an scv
     *
     * \note In the default case, we don't reconstruct anything. We do
     *       Reconstruction is only done when using the box method
     *       and enableReconstruction = true.
     *
     * \param spatialParams Class encapsulating the spatial parameters
     * \param element The finite element the scv is embedded in
     * \param scv The sub-control volume for which the saturation is computed
     * \param elemSol The solution at all dofs inside this element
     * \param sn The nonwetting phase saturation at the global dof
     */
    template<class SpatialParams, class Element, class Scv, class ElemSol>
    static typename ElemSol::PrimaryVariables::value_type
    reconstructSn(const SpatialParams& spatialParams,
                  const Element& element,
                  const Scv& scv,
                  const ElemSol& elemSol,
                  typename ElemSol::PrimaryVariables::value_type sn)
    { return sn; }
};

//! Specialization for the box scheme with the interface solver enabled
template<>
class TwoPScvSaturationReconstruction<DiscretizationMethod::box, /*enableReconstruction*/true>
{
public:
    /*!
     * \brief Compute the nonwetting phase saturation in an scv
     *
     * \param spatialParams Class encapsulating the spatial parameters
     * \param element The finite element the scv is embedded in
     * \param scv The sub-control volume for which the saturation is computed
     * \param elemSol The solution at all dofs inside this element
     * \param sn The nonwetting phase saturation at the global dof
     */
    template<class SpatialParams, class Element, class Scv, class ElemSol>
    static typename ElemSol::PrimaryVariables::value_type
    reconstructSn(const SpatialParams& spatialParams,
                  const Element& element,
                  const Scv& scv,
                  const ElemSol& elemSol,
                  typename ElemSol::PrimaryVariables::value_type sn)
    {
        // if this dof doesn't lie on a material interface, simply return Sn
        const auto& materialInterfaces = spatialParams.materialInterfaces();
        if (!materialInterfaces.isOnMaterialInterface(scv))
            return sn;

        // compute capillary pressure using material parameters associated with the dof
        const auto& interfacePcSw = materialInterfaces.pcSwAtDof(scv);
        const auto pc = interfacePcSw.pc(/*ww=*/1.0 - sn);

        // reconstruct by inverting the pc-sw curve
        const auto& pcSw = spatialParams.fluidMatrixInteraction(element, scv, elemSol).pcSwCurve();
        const auto pcMin = pcSw.endPointPc();

        if (pc < pcMin && pcMin > 0.0)
            return 0.0;
        else
            return 1.0 - pcSw.sw(pc);
    }
};

} // end namespace Dumux

#endif
