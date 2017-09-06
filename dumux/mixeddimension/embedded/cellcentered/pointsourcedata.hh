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
 * \ingroup EmbeddedCoupling
 * \brief Data associated with a point source
 */

#ifndef DUMUX_POINTSOURCEDATA_HH
#define DUMUX_POINTSOURCEDATA_HH

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(PrimaryVariables);
} // namespace Properties


/*!
 * \ingroup EmbeddedCoupling
 * \brief A point source data class used for integration in multidimension models
 * \note The point source and related data are connected via an identifier (id)
 */
template<class TypeTag>
class PointSourceData
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename std::vector<Dune::FieldVector<Scalar, 1> > ShapeValues;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables) BulkPrimaryVariables;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables) LowDimPrimaryVariables;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, SolutionVector) BulkSolutionVector;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, SolutionVector) LowDimSolutionVector;

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

public:
    void addBulkInterpolation(const ShapeValues& shapeValues,
                              const std::vector<unsigned int>& cornerIndices,
                              unsigned int eIdx)
    {
        bulkElementIdx_ = eIdx;
        bulkCornerIndices_ = cornerIndices;
        bulkShapeValues_ = shapeValues;
    }

    void addLowDimInterpolation(const ShapeValues& shapeValues,
                                const std::vector<unsigned int>& cornerIndices,
                                unsigned int eIdx)
    {
        lowDimElementIdx_ = eIdx;
        lowDimCornerIndices_ = cornerIndices;
        lowDimShapeValues_ = shapeValues;
    }

    void addBulkInterpolation(unsigned int eIdx)
    {
        bulkElementIdx_ = eIdx;
    }

    void addLowDimInterpolation(unsigned int eIdx)
    {
        lowDimElementIdx_ = eIdx;
    }

    BulkPrimaryVariables interpolateBulk(const BulkSolutionVector& sol)
    {
        BulkPrimaryVariables bulkPriVars(0.0);
        if (bulkIsBox)
        {
            for (unsigned int i = 0; i < bulkCornerIndices_.size(); ++i)
                for (unsigned int priVarIdx = 0; priVarIdx < bulkPriVars_.size(); ++priVarIdx)
                    bulkPriVars_[priVarIdx] += sol[bulkCornerIndices_[i]][priVarIdx]*bulkShapeValues_[i];
        }
        else
        {
            bulkPriVars = sol[bulkElementIdx()];
        }
        return bulkPriVars;
    }

    LowDimPrimaryVariables interpolateLowDim(const LowDimSolutionVector& sol)
    {
        LowDimPrimaryVariables lowDimPriVars(0.0);
        if (lowDimIsBox)
        {
            for (unsigned int i = 0; i < lowDimCornerIndices_.size(); ++i)
                for (unsigned int priVarIdx = 0; priVarIdx < lowDimPriVars_.size(); ++priVarIdx)
                    lowDimPriVars_[priVarIdx] += sol[lowDimCornerIndices_[i]][priVarIdx]*lowDimShapeValues_[i];
        }
        else
        {
            lowDimPriVars = sol[lowDimElementIdx()];
        }
        return lowDimPriVars;
    }

    unsigned int lowDimElementIdx() const
    { return lowDimElementIdx_; }

    unsigned int bulkElementIdx() const
    { return bulkElementIdx_; }

private:
    ShapeValues bulkShapeValues_, lowDimShapeValues_;
    std::vector<unsigned int> bulkCornerIndices_, lowDimCornerIndices_;
    BulkPrimaryVariables bulkPriVars_;
    LowDimPrimaryVariables lowDimPriVars_;
    unsigned int lowDimElementIdx_;
    unsigned int bulkElementIdx_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief A point source data class used for integration in multidimension models
 * \note The point source and related data are connected via an identifier (id)
 * When explicitly computing the circle average, i.e. the pressure for
 * the source term is computed as an integral over the circle describing
 * the surface of the one-dimensional tube. This exact determination of the bulk pressure
 * is necessary for convergence studies.
 */
template<class TypeTag>
class PointSourceDataCircleAverage : public PointSourceData<TypeTag>
{
    using ParentType = PointSourceData<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename std::vector<Dune::FieldVector<Scalar, 1> > ShapeValues;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, PrimaryVariables) BulkPrimaryVariables;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, PrimaryVariables) LowDimPrimaryVariables;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, SolutionVector) BulkSolutionVector;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, SolutionVector) LowDimSolutionVector;

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

public:
    PointSourceDataCircleAverage() : enableBulkCircleInterpolation_(false) {}

    BulkPrimaryVariables interpolateBulk(const BulkSolutionVector& sol)
    {
        // bulk interpolation on the circle is only enabled for source in the
        // lower dimensional domain if we use a circle distributed bulk sources
        // (see coupling manager circle sources)
        BulkPrimaryVariables bulkPriVars(sol[0]);
        bulkPriVars = 0.0;
        if (enableBulkCircleInterpolation_)
        {
            // compute the average of the bulk privars over the circle around the integration point
            // this computes $\bar{p} = \frac{1}{2\pi R} int_0^2*\pi p R \text{d}\theta.
            if (bulkIsBox)
            {
                assert(circleCornerIndices_.size() == circleShapeValues_.size());

                Scalar weightSum = 0.0;
                for (unsigned int j = 0; j < circleStencil_.size(); ++j)
                {
                    BulkPrimaryVariables priVars(0.0);
                    for (unsigned int i = 0; i < circleCornerIndices_[j].size(); ++i)
                        for (unsigned int priVarIdx = 0; priVarIdx < priVars.size(); ++priVarIdx)
                            priVars[priVarIdx] += sol[circleCornerIndices_[j][i]][priVarIdx]*circleShapeValues_[j][i];
                    // multiply with weight and add
                    priVars *= circleIpWeight_[j];
                    weightSum += circleIpWeight_[j];
                    bulkPriVars += priVars;
                }
                bulkPriVars /= weightSum;
            }
            else
            {
                Scalar weightSum = 0.0;
                for (unsigned int j = 0; j < circleStencil_.size(); ++j)
                {
                    for (unsigned int priVarIdx = 0; priVarIdx < bulkPriVars.size(); ++priVarIdx)
                        bulkPriVars[priVarIdx] += sol[circleStencil_[j]][priVarIdx]*circleIpWeight_[j];

                    weightSum += circleIpWeight_[j];
                }
                bulkPriVars /= weightSum;
            }
            return bulkPriVars;
        }
        else
        {
            return ParentType::interpolateBulk(sol);
        }
    }

    void addCircleInterpolation(const std::vector<std::vector<unsigned int> >& circleCornerIndices,
                                const std::vector<ShapeValues>& circleShapeValues,
                                const std::vector<Scalar>& circleIpWeight,
                                const std::vector<unsigned int>& circleStencil)
    {
        circleCornerIndices_ = circleCornerIndices;
        circleShapeValues_ = circleShapeValues;
        circleIpWeight_ = circleIpWeight;
        circleStencil_ = circleStencil;
        enableBulkCircleInterpolation_ = true;

    }

    void addCircleInterpolation(const std::vector<Scalar>& circleIpWeight,
                                const std::vector<unsigned int>& circleStencil)
    {
        circleIpWeight_ = circleIpWeight;
        circleStencil_ = circleStencil;
        enableBulkCircleInterpolation_ = true;
    }

    const std::vector<unsigned int>& circleStencil()
    {
        return circleStencil_;
    }

private:
    std::vector<std::vector<unsigned int> > circleCornerIndices_;
    std::vector<ShapeValues> circleShapeValues_;
    std::vector<Scalar> circleIpWeight_;
    std::vector<unsigned int> circleStencil_;
    bool enableBulkCircleInterpolation_;
};

} // end namespace Dumux

#endif
