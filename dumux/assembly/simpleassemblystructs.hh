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
 * \ingroup StaggeredDiscretization
 * \ingroup Assembly
 * \brief assembly structs
 */
#ifndef DUMUX_STAGGERED_SIMPLE_ASSEMBLY_STRUCTS_HH
#define DUMUX_STAGGERED_SIMPLE_ASSEMBLY_STRUCTS_HH

namespace Dumux {

/*!
* \ingroup StaggeredDiscretization
* \brief
*/
template<class CellCenterPrimaryVariables, class TypeTag>
struct SimpleMassBalanceSummands
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;

    CellCenterPrimaryVariables RHS;
    std::vector<CellCenterPrimaryVariables> coefficients;

    //! Constructor
    SimpleMassBalanceSummands(const Element& element, const FVElementGeometry& fvGeometry){
        this -> resize_(element, fvGeometry);
    }

    void setToZero(const Element& element, const FVElementGeometry& fvGeometry){
        this -> RHS = 0.0;

        for (auto&& scvf : scvfs(fvGeometry)){
            this -> coefficients[scvf.localFaceIdx()] = 0.0;
        }
    }

private:
    //I cannot simply resize a member without creating the object first
    void resize_(const Element& element, const FVElementGeometry& fvGeometry){
        const int numLocalFaces = element.subEntities(1);
        this -> coefficients.resize(numLocalFaces);
    }
};

/*!
* \ingroup StaggeredDiscretization
* \ingroup Assembly
* \brief
*/
template<class FacePrimaryVariables, class TypeTag>
struct SimpleMomentumBalanceSummands
{
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    FacePrimaryVariables RHS;
    FacePrimaryVariables selfCoefficient;
    FacePrimaryVariables oppositeCoefficient;
    FacePrimaryVariables pressureCoefficient;
    std::vector<FacePrimaryVariables> parallelCoefficients;
    std::vector<FacePrimaryVariables> innerNormalCoefficients;
    std::vector<FacePrimaryVariables> outerNormalCoefficients;

    //! Constructor
    SimpleMomentumBalanceSummands(const SubControlVolumeFace& scvf){
        this -> resize_(scvf);
    }

    void setToZero(const SubControlVolumeFace& scvf){
        const int numSubFaces = scvf.pairData().size();
        this -> RHS = 0.0;
        this -> selfCoefficient = 0.0;
        this -> oppositeCoefficient = 0.0;
        this -> pressureCoefficient = 0.0;
        for (int i = 0; i < numSubFaces; ++i){
            this -> parallelCoefficients[i] = 0.0;
            this -> innerNormalCoefficients[i] = 0.0;
            this -> outerNormalCoefficients[i] = 0.0;
        }
    }

private:
    void resize_(const SubControlVolumeFace& scvf){
        const int numSubFaces = scvf.pairData().size();
        this -> parallelCoefficients.resize(numSubFaces);
        this -> innerNormalCoefficients.resize(numSubFaces);
        this -> outerNormalCoefficients.resize(numSubFaces);
    }
};
} // namespace Dumux

#endif
