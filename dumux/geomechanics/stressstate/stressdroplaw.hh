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
 * \ingroup Geomechanics
 * \brief Stress drop law checks if the shear failure occurs and reduces
 *  the shear stress directly from the failure stress state.
 *  The failure stress state maybe rotated.
 */
#ifndef DUMUX_STRESS_STATE_STRESS_DROP_LAW_HH
#define DUMUX_STRESS_STATE_STRESS_DROP_LAW_HH

#include <cmath>
#include <algorithm>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dumux/common/math.hh>
#include "math/mohrspace.hh"

namespace Dumux
{

    template <class Scalar, class StressTensor, class MohrSpaceTypeTraits>
    class StressDropLaw
    {
        using Point = typename MohrSpaceTypeTraits::PointType;
        using MohrCircle = typename MohrSpaceTypeTraits::MohrCircleType;


        static constexpr int dim = StressTensor::rows;
        static_assert(dim == 2, "Stress drop law works only for the 2D tensor.");

    public:
        StressDropLaw(const StressTensor &stress)
            : stress_(stress)
            , mohrCircle_(stress)
            , stressPoint_(getStressPoint_(stress))
            , isInUpperHalfPlane_(stressPoint_.y() > 0.0 ? true : false)
        {}

        template <class StressLawParams>
        bool hasFailure( const StressLawParams& params) const
        { return mohrCircle_.hasIntersection(params.envelopeLine(isInUpperHalfPlane_));}

        template <class StressDropLawParams>
        void doStressDrop(const StressDropLawParams &param)
        {
            const auto envelopeLine = param.envelopeLine(isInUpperHalfPlane_);
            auto intersections = mohrCircle_.intersections(envelopeLine);

            // x0 >= x1, and within the range, the stress drop can be done directly
            if (stressPoint_.x() == std::clamp(stressPoint_.x(), intersections[1].x(), intersections[0].x()))
            {
                doStressDrop_(param);
            }

            // x < x1 or x > x0
            else
            {
                const auto closetIntersection = stressPoint_.x() > intersections[0].x() ? intersections[0] : intersections[1];

                /// calculate angle
                Scalar angle = angleBtwVectors(closetIntersection,
                stressPoint_, mohrCircle_.center());

                if (!isInUpperHalfPlane_){ angle*= -1.0;}
                if (stressPoint_.x() < intersections[0].x()){angle*= -1.0;}

                /// rotate in the physical space
                rotateStressTensor(angle/2);

                ///
                doStressDrop_(param);

                /// rotate back to the orthogonal plane
                rotateStressTensor(-angle/2);
            }
        }

        /*!
         * \brief calculate the angle between two vectors
         * OP1 and OP2
         *
         * \param p1 Point p1
         * \param p2 Point p2
         * \param pO Point pO
         * \return Scalar angle in radian
         * \note angle in Mohr space, \f$ 2 \theta\f$
         */
        Scalar angleBtwVectors(const Point& p1,
                               const Point& p2,
                               const Point& pO) const
        {
            const auto op1 = p1.pos() - pO.pos();
            const auto op2 = p2.pos() - pO.pos();
            using std::acos;
            Scalar angle = acos(op1 * op2 / op1.two_norm() / op2.two_norm());
            return angle;
        }

        StressTensor stressTensor() const
        {
            return stress_;
        }

    private:
        Point getStressPoint_(const StressTensor& stress)
        {
            if (stress[0][0] > stress[1][1])
            {
                return Point(stress[0][0], -stress[0][1]);
            }
            else
            {
                return Point(stress[1][1], stress[0][1]);
            }
        }

        /*!
         * \brief calculate the transformation of stress tensor
         *        Convention: counter-clockwise rotation is positive
         *        see also Figure 5 Signconvention 3 in
         *        https://en.wikipedia.org/wiki/Mohr%27s_circle.
         *
         * \param angle angle in radian
         * \note angle in Mohr Space, \f$ 2 \theta \f$
         */
        void rotateStressTensor(const Scalar &angle)
        {
            StressTensor rotationTensor(0.0);
            using std::cos;
            using std::sin;
            rotationTensor[0][0] = cos(angle);
            rotationTensor[0][1] = -sin(angle);
            rotationTensor[1][0] = sin(angle);
            rotationTensor[1][1] = cos(angle);

            stress_ = stress_.leftmultiplyany(getTransposed(rotationTensor));
            stress_ = stress_.rightmultiplyany(rotationTensor);
        }

        /*!
         * \brief reduce the shear stress in stress tensor.
         *
         * \param params
         */
        template <class StressDropLawParams>
        void doStressDrop_(const StressDropLawParams &params)
        {
            if (std::abs(stress_[0][1]) < params.stressDrop())
            {
                stress_[0][1] = stress_[1][0] = 0.0;
                return;
            }

            // the absolute value of shear stress should be smaller
            if (stress_[0][1] > 0)
            {
                stress_[0][1] -= params.stressDrop();
                stress_[1][0] -= params.stressDrop();
            }
            else
            {
                stress_[0][1] += params.stressDrop();
                stress_[1][0] += params.stressDrop();
            }
        }

        StressTensor stress_;
        MohrCircle mohrCircle_;
        Point stressPoint_; // the concerned stress point

        const bool isInUpperHalfPlane_;
    };// end class StressDropLaw

} // namespace Dumux
#endif
