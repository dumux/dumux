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
 * \brief A newton solver specific to the Richards problem.
 */
#ifndef DUMUX_RICHARDS_NEWTON_CONTROLLER_HH
#define DUMUX_RICHARDS_NEWTON_CONTROLLER_HH

#include "properties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup Newton
 * \brief A Richards model specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton controller.
 */
template <class TypeTag>
class RichardsNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pwIdx = Indices::pwIdx };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
    static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead);

public:
    /*!
     * \brief Constructor
     */
    RichardsNewtonController(const Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Update the current solution of the newton method
     *
     * This is basically the step
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        ParentType::newtonUpdate(uCurrentIter, uLastIter, deltaU);

        if ( (!GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, UseLineSearch)) && (!useHead) )
        {
            // do not clamp anything after 5 iterations
            if (this->numSteps_ > 4)
                return;

            // clamp saturation change to at most 20% per iteration
            FVElementGeometry fvGeometry;
            const GridView &gridView = this->problem_().gridView();
            for (const auto& element : elements(gridView)) {
                fvGeometry.update(gridView, element);
                using std::max;
                using std::min;
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int dofIdxGlobal = this->model_().dofMapper().subIndex(element, scvIdx, dofCodim);

                    // calculate the old wetting phase saturation
                    const SpatialParams &spatialParams = this->problem_().spatialParams();
                    const MaterialLawParams &mp = spatialParams.materialLawParams(element, fvGeometry, scvIdx);
                    Scalar pcMin = MaterialLaw::pc(mp, 1.0);
                    Scalar pw = uLastIter[dofIdxGlobal][pwIdx];
                    Scalar pn = max(this->problem_().referencePressure(element, fvGeometry, scvIdx),
                                         pw + pcMin);
                    Scalar pcOld = pn - pw;
                    Scalar SwOld = max<Scalar>(0.0, MaterialLaw::sw(mp, pcOld));

                    // convert into minimum and maximum wetting phase
                    // pressures
                    Scalar pwMin = pn - MaterialLaw::pc(mp, SwOld - 0.2);
                    Scalar pwMax = pn - MaterialLaw::pc(mp, SwOld + 0.2);

                    // clamp the result
                    pw = uCurrentIter[dofIdxGlobal][pwIdx];
                    pw = max(pwMin, min(pw, pwMax));
                    uCurrentIter[dofIdxGlobal][pwIdx] = pw;

                }
            }
        }
    }
};
}

#endif
