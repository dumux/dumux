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
// class DecoupledTwoPElementVolumeVariables : public CCElementVolumeVariables<TypeTag>
class DecoupledTwoPElementVolumeVariables : public BoxElementVolumeVariables<TypeTag>
{
//     typedef CCElementVolumeVariables<TypeTag> ParentType;
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

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief The constructor.
     */
    DecoupledTwoPElementVolumeVariables()
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

//     /*!
//      * \brief Construct the volume variables for all of vertices of an element.
//      *
//      * \param problem The problem which needs to be simulated.
//      * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
//      * \param fvGeometry The finite volume geometry of the element
//      * \param oldSol Tells whether the model's previous or current solution should be used.
//      */
//
//     void update(const SolutionVector &globalSol,
//                 const Problem &problem,
//                 const Element &element,
//                 const FVElementGeometry &fvGeometry,
//                 const bool oldSol)
//     {
//         const VertexMapper &vertexMapper = problem.vertexMapper();
//         // we assert that the i-th shape function is
//         // associated to the i-th vertex of the element.
//         int numVertices = element.subEntities(dim);
//         this->resize(numVertices);
//         for (int scvIdx = 0; scvIdx < numVertices; scvIdx++) {
//             const PrimaryVariables &priVars
//                 = globalSol[vertexMapper.subIndex(element, scvIdx, dim)];
//
//             // reset evaluation point to zero
//             (*this)[scvIdx].setEvalPoint(0);
//
//             (*this)[scvIdx].update(priVars,
//                               problem,
//                               element,
//                               fvGeometry,
//                               scvIdx,
//                               oldSol);
//         }
//
//         this->updateEffPorosity(problem, element, fvGeometry, oldSol);
//     }


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
        int numScv = fvGeometry.numScv;

//         std::cout << "elemVolVars: numScv = " << numScv << std::endl;

        typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
        typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        (*this)[scvIdx].effPorosity_ = (*this)[scvIdx].initialPorosity();

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        {

                LocalPosition scvCenter = fvGeometry.subContVol[scvIdx].localCenter;
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


            for (int j = 0; j < shapeVal.size(); ++j)
            {
// //                 B += problem.getBNodeWiseAveraged(element, fvGeometry, j) * shapeVal[j];
// //                B += shapeVal[j]/problem.getBNodeWiseAveraged(element, fvGeometry, j);
                pWCurrentIteration +=  (*this)[j].pressure(wPhaseIdx) * shapeVal[j];
                pNCurrentIteration +=  (*this)[j].pressure(nPhaseIdx) * shapeVal[j];
                sWCurrentIteration +=  (*this)[j].saturation(wPhaseIdx) * shapeVal[j];
                sNCurrentIteration +=  (*this)[j].saturation(nPhaseIdx) * shapeVal[j];

                // for fixed stress, getpWOldIteration_() in the problem  returns p and S of old iteration
                // for pore compressibility getpWOldIteration_() in the problem  returns pInit and Sinit
                pWOldIteration +=  problem.getpWOldIteration_(element, fvGeometry, j) * shapeVal[j];
                pNOldIteration +=  problem.getpNOldIteration_(element, fvGeometry, j) * shapeVal[j];
                sWOldIteration +=  problem.getsWOldIteration_(element, fvGeometry, j) * shapeVal[j];
                sNOldIteration +=  problem.getsNOldIteration_(element, fvGeometry, j) * shapeVal[j];
            }


//                 pWCurrentIteration =  (*this)[scvIdx].pressure(wPhaseIdx);
//                 pNCurrentIteration =  (*this)[scvIdx].pressure(nPhaseIdx);
//                 sWCurrentIteration =  (*this)[scvIdx].saturation(wPhaseIdx);
//                 sNCurrentIteration =  (*this)[scvIdx].saturation(nPhaseIdx);
//
//                 pWOldIteration =  problem.getpWOldIteration_(element, fvGeometry, scvIdx);
//                 pNOldIteration =  problem.getpNOldIteration_(element, fvGeometry, scvIdx);
//                 sWOldIteration =  problem.getsWOldIteration_(element, fvGeometry, scvIdx);
//                 sNOldIteration =  problem.getsNOldIteration_(element, fvGeometry, scvIdx);

            (*this)[scvIdx].deltaVolumetricStrainOldIteration_ = problem.getDeltaVolumetricStrainOldIteration(element, fvGeometry, scvIdx);

            Scalar effPressureCurrentIteration = pWCurrentIteration *  sWCurrentIteration +
                                                    pNCurrentIteration *  sNCurrentIteration;

            Scalar effPressureOldIteration = pWOldIteration *  sWOldIteration +
                                                    pNOldIteration *  sNOldIteration;

            // for fixed-stress
            (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_ = 1.0/B * (effPressureCurrentIteration -   effPressureOldIteration) + (*this)[scvIdx].deltaVolumetricStrainOldIteration_;
            // for pore compressibility
//             (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_ = (*this)[scvIdx].initialPorosity() * poreCompressibility * (effPressureCurrentIteration -   effPressureOldIteration);

            if(problem.coupled() == true)
            {
                if (isOldSol == true)
                {
                    (*this)[scvIdx].effPorosity_ =  problem.getEffPorosityOldTimestep(element, fvGeometry, scvIdx);
                }
                else
                {
                    if(problem.evalOriginalRhs() == true)
                    {
                        Scalar dampingFactor = GET_RUNTIME_PARAM(TypeTag, Scalar,Damping.DampingFactor);

                        (*this)[scvIdx].effPorosity_ =
                        dampingFactor * ((*this)[scvIdx].initialPorosity() + (*this)[scvIdx].deltaVolumetricStrainOldIteration_)
                        + (1.0 - dampingFactor) * ((*this)[scvIdx].initialPorosity() + (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_);
//                         std::cout << "deltaVolumetricStrainOldIteration_[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].deltaVolumetricStrainOldIteration_ << std::endl;
                    }
                    if(problem.evalGlobalResidual() == true)
                    {
                        (*this)[scvIdx].effPorosity_ =
                        ((*this)[scvIdx].initialPorosity() + (*this)[scvIdx].deltaVolumetricStrainOldIteration_);
                    }
                    else
                    {
                        //                     (*this)[scvIdx].effPorosity_ = 1 - (1 - (*this)[scvIdx].initialPorosity() )*exp( -((*this)[scvIdx].deltaVolumetricStrainCurrentIteration_));
                        (*this)[scvIdx].effPorosity_ = ((*this)[scvIdx].initialPorosity() + (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_)/*/(1.0 + (*this)[scvIdx].deltaVolumetricStrainCurrentIteration_)*/;
/*
                        std::cout << "effPorosityNew_[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].effPorosity_ << std::endl;*/
                    }
                }
            }
        }
    }
};

} // namespace Dumux

#endif
