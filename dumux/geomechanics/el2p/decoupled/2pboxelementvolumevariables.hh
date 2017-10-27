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
 * \brief Volume variables gathered on an element
 */
#ifndef DUMUX_DECOUPLED_TWOP_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_DECOUPLED_TWOP_ELEMENT_VOLUME_VARIABLES_HH

#include <dumux/implicit/box/elementvolumevariables.hh>


namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class DecoupledTwoPBoxElementVolumeVariables : public BoxElementVolumeVariables<TypeTag>
{
    typedef BoxElementVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

public:
    /*!
     * \brief The constructor.
     */
    DecoupledTwoPBoxElementVolumeVariables()
    { }

    /*!
     * \brief Construct the volume variables for all of vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param oldSol Tells whether the model's previous or current solution should be used.
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const bool oldSol)
    {
        ParentType::update(problem,
                           element,
                           fvGeometry,
                           oldSol);

        this->updateEffPorosity(problem, element, fvGeometry, oldSol);
    }

    /*!
     * \brief Construct the volume variables for all of vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param oldSol Tells whether the model's previous or current solution should be used.
     */

    void update(const SolutionVector &globalSol,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const bool oldSol)
    {
        const VertexMapper &vertexMapper = problem.vertexMapper();
        // we assert that the i-th shape function is
        // associated to the i-th vertex of the element.
        int numVertices = element.subEntities(dim);
        this->resize(numVertices);
        for (int scvIdx = 0; scvIdx < numVertices; scvIdx++) {
            const PrimaryVariables &priVars
                = globalSol[vertexMapper.subIndex(element, scvIdx, dim)];

            // reset evaluation point to zero
            (*this)[scvIdx].setEvalPoint(0);

            (*this)[scvIdx].update(priVars,
                              problem,
                              element,
                              fvGeometry,
                              scvIdx,
                              oldSol);
        }

        this->updateEffPorosity(problem, element, fvGeometry, oldSol);
    }


    /*!
     * \brief Construct the volume variables for all of vertices of an
     *        element given a solution vector computed by PDELab.
     *
     * \tparam ElementSolutionVector The container type which stores the
     *                           primary variables of the element
     *                           using _local_ indices
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param elementSolVector The local solution for the element using PDELab ordering
     */
    template<typename ElementSolutionVector>
    void updatePDELab(const Problem &problem,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const ElementSolutionVector& elementSolVector)
    {
        ParentType::update(problem,
                           element,
                           fvGeometry,
                           elementSolVector);
    }

    /*!
     * \brief Update the effective porosities for all vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param isOldSol Specifies whether this is the previous solution or the current one
     *
     * This function is required for the update of the effective porosity values at the
     * vertices.
     *
     * During the partial derivative calculation, changes of the solid displacement
     * at vertex i can affect effective porosities of all element vertices.
     * To correctly update the effective porosities of all element vertices
     * an iteration over all scv faces is required.
     * The remaining volvars are only updated for the vertex whose primary variable
     * is changed for the derivative calculation.
     */
    void updateEffPorosity(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                bool isOldSol)
    {
        int eIdx = problem.model().elementMapper().index(element);

        int numScv = element.subEntities(dim);

        typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
        typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        (*this)[scvIdx].effPorosity_ = (*this)[scvIdx].initialPorosity();

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        {
            LocalPosition scvCenter = fvGeometry.subContVol[scvIdx].localCenter;
            GlobalPosition scvCenterGlobal = element.geometry().global(fvGeometry.subContVol[scvIdx].localCenter);
//             std::cout << "scvCenterGlobal = " << scvCenterGlobal << std::endl;
            const LocalFiniteElementCache feCache;
            Dune::GeometryType geomType = element.geometry().type();

            const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
            std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
                        localFiniteElement.localBasis().evaluateFunction(scvCenter, shapeVal);
            Scalar pWCurrentIteration =  0.0;
            Scalar pNCurrentIteration =  0.0;
            Scalar sWCurrentIteration =  0.0;
            Scalar sNCurrentIteration =  0.0;

            Scalar pWOldIteration =  0.0;
            Scalar pNOldIteration =  0.0;
            Scalar sWOldIteration =  0.0;
            Scalar sNOldIteration =  0.0;

            // drained bulk modulus
            Scalar B = 0.0;
            B = problem.spatialParams().lameParams(element, fvGeometry, scvIdx)[1];
            Scalar poreCompressibility = problem.spatialParams().lameParams(element, fvGeometry, scvIdx)[3];
//                B = problem.getBNodeWiseAveraged(element, fvGeometry, scvIdx);

//             if (eIdx == 211)
//                 std::cout << "pw[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].pressure(wPhaseIdx) << std::endl;

//             std::cout << "problem.getpWOldIteration_(element, fvGeometry," << scvIdx << ") at element " << eIdx << " = " << problem.getpWOldIteration_(element, fvGeometry, scvIdx) << std::endl;
//             std::cout << "(*this)[" << scvIdx << "].pressure(wPhaseIdx) = " << (*this)[scvIdx].pressure(wPhaseIdx) << std::endl;

            for (int j = 0; j < shapeVal.size(); ++j)
            {
//                 B += problem.getBNodeWiseAveraged(element, fvGeometry, j) * shapeVal[j];
//                B += shapeVal[j]/problem.getBNodeWiseAveraged(element, fvGeometry, j);
                pWCurrentIteration +=  (*this)[j].pressure(wPhaseIdx) * shapeVal[j];
                pNCurrentIteration +=  (*this)[j].pressure(nPhaseIdx) * shapeVal[j];
                sWCurrentIteration +=  (*this)[j].saturation(wPhaseIdx) * shapeVal[j];
                sNCurrentIteration +=  (*this)[j].saturation(nPhaseIdx) * shapeVal[j];

                // for fixed stress
                pWOldIteration +=  problem.getpWOldIteration_(element, fvGeometry, j) * shapeVal[j];
                pNOldIteration +=  problem.getpNOldIteration_(element, fvGeometry, j) * shapeVal[j];
                sWOldIteration +=  problem.getsWOldIteration_(element, fvGeometry, j) * shapeVal[j];
                sNOldIteration +=  problem.getsNOldIteration_(element, fvGeometry, j) * shapeVal[j];

//                 // for pore compressibility
//                 pWOldIteration +=  problem.getpInitPerScv_(element, fvGeometry, j) * shapeVal[j];
//                 sWOldIteration =  1.0;


//                 pWCurrentIteration +=  (*this)[scvIdx].pressure(wPhaseIdx);
//                 pNCurrentIteration +=  (*this)[scvIdx].pressure(nPhaseIdx);
//                 sWCurrentIteration +=  (*this)[scvIdx].saturation(wPhaseIdx);
//                 sNCurrentIteration +=  (*this)[scvIdx].saturation(nPhaseIdx);
//
//                 pWOldIteration +=  problem.getpWOldIteration_(element, fvGeometry, scvIdx);
//                 pNOldIteration +=  problem.getpNOldIteration_(element, fvGeometry, scvIdx);
//                 sWOldIteration +=  problem.getsWOldIteration_(element, fvGeometry, scvIdx);
//                 sNOldIteration +=  problem.getsNOldIteration_(element, fvGeometry, scvIdx);
            }


//             B = 1.0 / B;

//                 if (eIdx == 0)
//                 {
//                     std::cout << "scvCenter_[" << eIdx << "][" << scvIdx << "] = " << scvCenter[0] << " " << scvCenter[1] << std::endl;
//
//                     for (int j = 0; j < shapeVal.size(); ++j)
//                     {
//                         std::cout << "elemVolVarsTranspProblem[" << j << "] = " << (*this)[j].pressure(wPhaseIdx) << std::endl;
//                         std::cout << "shapeVal[" << j << "] = " << shapeVal[j] << std::endl;
//                     }
//                 }

//             Scalar deltaVolumetricStrainOldIteration = 0.0;

            (*this)[scvIdx].deltaVolumetricStrainOldIteration_ = problem.getDeltaVolumetricStrainOldIteration(element, fvGeometry, scvIdx);

//             Scalar eps = 1e-6;
//             if(problem.coupled() == true)
//             if ((scvCenterGlobal[0] > 500 - eps) && (scvCenterGlobal[0] < 1500 + eps) &&
//                 (scvCenterGlobal[1] > 500 - eps) && (scvCenterGlobal[1] < 1500 + eps)
//             )
//             {
//                 std::cout.precision(15);
//                 std::cout << "pWCurrentIteration[" << eIdx << "][" << scvIdx << "] = " << pWCurrentIteration << std::endl;
//                 std::cout << "deltaVolumetricStrainOldIteration[" << eIdx << "][" << scvIdx << "] = " << deltaVolumetricStrainOldIteration << std::endl;
//             }

            Scalar effPressureCurrentIteration = pWCurrentIteration *  sWCurrentIteration +
                                                    pNCurrentIteration *  sNCurrentIteration;

            Scalar effPressureOldIteration = pWOldIteration *  sWOldIteration +
                                                    pNOldIteration *  sNOldIteration;


            Scalar pressureDifference = effPressureCurrentIteration - effPressureOldIteration;


//             if( eIdx == 12873)
//             {
//             if (std::abs(pressureDifference/effPressureCurrentIteration) > 1.0e-5)
//             {
//                 std::cout.precision(15);
//                 std::cout << "pWCurrentIteration[" << eIdx << "][" << scvIdx << "] = " << pWCurrentIteration << std::endl;
//                 std::cout << "effPressureCurrentIteration[" << eIdx << "][" << scvIdx << "] = " << effPressureCurrentIteration << std::endl;
//                 std::cout << "effPressureOldIteration[" << eIdx << "][" << scvIdx << "] = " << effPressureOldIteration << std::endl;
//
//                 std::cout << "B[" << eIdx << "][" << scvIdx << "] = " << B << std::endl;

//                 std::cout << "deltaVolumetricStrainOldIteration[" << eIdx << "][" << scvIdx << "] = " << deltaVolumetricStrainOldIteration << std::endl;
//             }
//             }

            // for fixed-stress
//             (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_ = 1.0/B * (effPressureCurrentIteration -   effPressureOldIteration) + (*this)[scvIdx].deltaVolumetricStrainOldIteration_;
            // for pore compressibility
//             (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_ = poreCompressibility * (effPressureCurrentIteration -   effPressureOldIteration);

//             if( eIdx == 211)
//             {
//                 std::cout << "pWOldIteration    [" << eIdx << "][" << scvIdx << "] = " << pWOldIteration << std::endl;
//                 std::cout << "pWCurrentIteration[" << eIdx << "][" << scvIdx << "] = " << pWCurrentIteration << std::endl;
//             }
//
//             std::cout << "deltaVolumetricStrainCurrentIteration = " << deltaVolumetricStrainCurrentIteration << std::endl;
//             std::cout << "deltaVolumetricStrainOldIteration = " << deltaVolumetricStrainOldIteration << std::endl;

            if(problem.coupled() == true)
            {
                if (isOldSol == true)
                {
//                     Scalar deltaVolumetricStrainOldTimestep = 1.0/B * (effPressureOldIteration -   effPressureOldIteration) + deltaVolumetricStrainOldIteration;
//                     (*this)[scvIdx].effPorosity_ = 1 - (1 - (*this)[scvIdx].initialPorosity() )*exp( -(deltaVolumetricStrainOldTimestep));
                    (*this)[scvIdx].effPorosity_ =  problem.getEffPorosityOldTimestep(element, fvGeometry, scvIdx);

//                     if( ((scvCenterGlobal[0] > -0.5) && (scvCenterGlobal[0] < 50)) &&
//                     ((scvCenterGlobal[1] > -0.1) && (scvCenterGlobal[1] < 50)) )
//                     {
//                         std::cout << "effPorosityOld_[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].effPorosity_ << std::endl;
//                     }

                }
                else
                {
                    (*this)[scvIdx].effPorosity_ = (*this)[scvIdx].initialPorosity() * (1.0 + poreCompressibility * pressureDifference);
//                     (*this)[scvIdx].effPorosity_ = problem.getEffPorosityOldTimestep(element, fvGeometry, scvIdx)/(1.0 - poreCompressibility * pressureDifference);
                    //                     (*this)[scvIdx].effPorosity_ = 1 - (1 - (*this)[scvIdx].initialPorosity() )*exp( -((*this)[scvIdx].deltaVolumetricStrainCurrentIteration_));
//                     (*this)[scvIdx].effPorosity_ = ((*this)[scvIdx].initialPorosity() + deltaVolumetricStrainCurrentIteration)/(1.0 + deltaVolumetricStrainCurrentIteration);

                    if( (eIdx == 211) && (scvIdx == 0) )
//                     if( ((scvCenterGlobal[0] > -0.5) && (scvCenterGlobal[0] < 50)) &&
//                     ((scvCenterGlobal[1] > -0.1) && (scvCenterGlobal[1] < 50)) )
                    {
//                         std::cout << "pWCurrentIteration[" << eIdx << "][" << scvIdx << "] = " << pWCurrentIteration << std::endl;
//                         std::cout << "pWOldIteration[" << eIdx << "][" << scvIdx << "] = " << pWOldIteration << std::endl;
//
// //                         std::cout << "deltaVolumetricStrainOldIteration[" << eIdx << "][" << scvIdx << "] = " << deltaVolumetricStrainOldIteration << std::endl;
//                         std::cout << "deltaVolumetricStrainCurrentIteration[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_ << std::endl;
//
//                         std::cout << "effPorosity_[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].effPorosity_ << std::endl;
                    }
                }
            }
        }
    }
};

} // namespace Dumux

#endif
