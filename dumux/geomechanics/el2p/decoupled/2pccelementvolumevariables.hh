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
#ifndef DUMUX_DECOUPLED_TWOP_CC_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_DECOUPLED_TWOP_CC_ELEMENT_VOLUME_VARIABLES_HH

#include <dumux/implicit/box/elementvolumevariables.hh>
#include <dumux/implicit/cellcentered/elementvolumevariables.hh>


namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class DecoupledTwoPCCElementVolumeVariables : public CCElementVolumeVariables<TypeTag>
{
    typedef CCElementVolumeVariables<TypeTag> ParentType;

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
    DecoupledTwoPCCElementVolumeVariables()
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

        (*this)[0].effPorosity_ = (*this)[0].initialPorosity();

        Scalar pWCurrentIteration =  0.0;
        Scalar pNCurrentIteration =  0.0;
        Scalar sWCurrentIteration =  0.0;
        Scalar sNCurrentIteration =  0.0;

        Scalar pWOldIteration =  0.0;
        Scalar pNOldIteration =  0.0;
        Scalar sWOldIteration =  0.0;
        Scalar sNOldIteration =  0.0;

        // drained bulk modulus
        // drained bulk modulus
        Scalar B = 0.0;
        Scalar lambda = problem.spatialParams().lameParams(element, fvGeometry, 0)[0];
        Scalar mu = problem.spatialParams().lameParams(element, fvGeometry, 0)[1];

        B = lambda + 2.0/3.0*mu;

        Scalar poreCompressibility = problem.spatialParams().lameParams(element, fvGeometry, 0)[2];

        pWCurrentIteration +=  (*this)[0].pressure(wPhaseIdx);
        pNCurrentIteration +=  (*this)[0].pressure(nPhaseIdx);
        sWCurrentIteration +=  (*this)[0].saturation(wPhaseIdx);
        sNCurrentIteration +=  (*this)[0].saturation(nPhaseIdx);

        pWOldIteration +=  problem.getpWOldIteration_(element, fvGeometry, 0);
        pNOldIteration +=  problem.getpNOldIteration_(element, fvGeometry, 0);
        sWOldIteration +=  problem.getsWOldIteration_(element, fvGeometry, 0);
        sNOldIteration +=  problem.getsNOldIteration_(element, fvGeometry, 0);

        (*this)[0].deltaVolumetricStrainOldIteration_ = problem.getDeltaVolumetricStrainOldIteration(element, fvGeometry, 0);

        Scalar effPressureCurrentIteration = pWCurrentIteration *  sWCurrentIteration +
                                                pNCurrentIteration *  sNCurrentIteration;

        Scalar effPressureOldIteration = pWOldIteration *  sWOldIteration +
                                                pNOldIteration *  sNOldIteration;

        // for pore compressibility
        if (GET_RUNTIME_PARAM(TypeTag, bool,PoreCompressibility.UsePoreCompressibility))
            (*this)[0].deltaVolumetricStrainCurrentIteration_ = (*this)[0].initialPorosity() * poreCompressibility * (effPressureCurrentIteration - effPressureOldIteration);
        // for fixed-stress
        else
            (*this)[0].deltaVolumetricStrainCurrentIteration_ = 1.0/B * (effPressureCurrentIteration -   effPressureOldIteration) + (*this)[0].deltaVolumetricStrainOldIteration_;

        if(problem.coupled() == true)
        {
            if (isOldSol == true)
            {
                (*this)[0].effPorosity_ =  problem.getEffPorosityOldTimestep(element, fvGeometry, 0);
            }
            else
            {
//                     //                     (*this)[0].effPorosity_ = 1 - (1 - (*this)[0].initialPorosity() )*exp( -((*this)[0].deltaVolumetricStrainCurrentIteration_));
                (*this)[0].effPorosity_ = ((*this)[0].initialPorosity() + (*this)[0].deltaVolumetricStrainCurrentIteration_)/*/(1.0 + (*this)[0].deltaVolumetricStrainCurrentIteration_)*/;
//
//                     if (eIdx == 2925)
//                     {
//                         std::cout << "effPressureCurrentIteration[" << eIdx << "][" << 0 << "] = " << effPressureCurrentIteration << std::endl;
//
//                         std::cout << "effPressureOldIteration[" << eIdx << "][" << 0 << "] = " << effPressureOldIteration << std::endl;
//
// //                         std::cout << "effPressureInit[" << eIdx << "][" << 0 << "] = " << effPressureInit << std::endl;
//
//                         std::cout << "deltaVolumetricStrainOldIteration_[" << eIdx << "][" << 0 << "] = " << (*this)[0].deltaVolumetricStrainOldIteration_ << std::endl;
//
//                         std::cout << "deltaVolumetricStrainCurrentIteration_[" << eIdx << "][" << 0 << "] = " << (*this)[0].deltaVolumetricStrainCurrentIteration_ << std::endl;
//                         std::cout << "effPorosityNew_[" << eIdx << "][" << 0 << "] = " << (*this)[0].effPorosity_ << std::endl;
//                     }
            }
        }

    }
};

} // namespace Dumux

#endif
