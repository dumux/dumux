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
#ifndef DUMUX_DECOUPLED_ELASTIC_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_DECOUPLED_ELASTIC_ELEMENT_VOLUME_VARIABLES_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <dumux/implicit/box/elementvolumevariables.hh>
#include "dumux/geomechanics/el2p/properties.hh"

namespace Dumux
{

/*!
 * \ingroup ElTwoPBoxModel
 *
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class DecoupledElasticElementVolumeVariables : public std::vector<typename GET_PROP_TYPE(TypeTag, VolumeVariables) >
{
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename Element::Geometry::JacobianInverseTransposed JacobianInverseTransposed;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;

    typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

public:
    /*!
     * \brief The constructor.
     */
    DecoupledElasticElementVolumeVariables()
    { }

    /*!
     * \brief Construct the volume variables for all vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param isOldSol Tells whether the model's previous or current solution should be used.
     *
     * This class is required for the update of the effective porosity values at the
     * vertices since it is a function of the divergence of the solid displacement
     * at the integration points
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                bool isOldSol)
    {
        // retrieve the current or the previous solution vector and write the values into globalSol
        const SolutionVector &globalSol =
            isOldSol?
            problem.model().prevSol():
            problem.model().curSol();

        const GridFunctionSpace& gridFunctionSpace = problem.model().jacobianAssembler().gridFunctionSpace();
        const typename GridFunctionSpace::Ordering& ordering = gridFunctionSpace.ordering();
        // copy the values of the globalSol vector to the localFunctionSpace values of the current element
        LocalFunctionSpace localFunctionSpace(gridFunctionSpace);
        localFunctionSpace.bind(element);
        std::vector<Scalar> values(localFunctionSpace.size());
        for (typename LocalFunctionSpace::Traits::IndexContainer::size_type k=0; k<localFunctionSpace.size(); ++k)
        {
            const typename GridFunctionSpace::Ordering::Traits::DOFIndex& di = localFunctionSpace.dofIndex(k);
            typename GridFunctionSpace::Ordering::Traits::ContainerIndex ci;
            ordering.mapIndex(di.view(),ci);
            values[k] = globalSol[ci];
        }

        typedef typename LocalFunctionSpace::template Child<0>::Type ScalarDispLFS;

        int numScv = element.subEntities(dim);
        this->resize(numScv);

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        {
            // solution vector solI for each vertex
            PrimaryVariables solI;
            // solid displacement values for each coordinate direction
            for (int coordDir = 0; coordDir < dim; coordDir++)
            {
                const ScalarDispLFS& scalarDispLFS = localFunctionSpace.child(coordDir);
                solI[Indices::u(coordDir)] = values[scalarDispLFS.localIndex(scvIdx)];
            }
            // reset evaluation point to zero
            (*this)[scvIdx].setEvalPoint(0);

            (*this)[scvIdx].update(solI,
                              problem,
                              element,
                              fvGeometry,
                              scvIdx,
                              isOldSol);

            Valgrind::CheckDefined((*this)[scvIdx]);
        }
        this->updateEffPorosity(problem, element, fvGeometry, isOldSol);

        if (isOldSol)
            prevValues_ = values;
        else
            dofValues_ = values;
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
        int numScv = element.subEntities(dim);

        // retrieve the current or the previous solution vector and write the values into globalSol
        const SolutionVector &globalSol =
            isOldSol?
            problem.model().prevSol():
            problem.model().curSol();

        // copy the values of the globalSol vector to the localFunctionSpace values of the current element
        const GridFunctionSpace& gridFunctionSpace = problem.model().jacobianAssembler().gridFunctionSpace();
        const typename GridFunctionSpace::Ordering& ordering = gridFunctionSpace.ordering();
        LocalFunctionSpace localFunctionSpace(gridFunctionSpace);
        localFunctionSpace.bind(element);
        std::vector<Scalar> values(localFunctionSpace.size());
        for (typename LocalFunctionSpace::Traits::IndexContainer::size_type k=0; k<localFunctionSpace.size(); ++k)
        {
            const typename GridFunctionSpace::Ordering::Traits::DOFIndex& di = localFunctionSpace.dofIndex(k);
            typename GridFunctionSpace::Ordering::Traits::ContainerIndex ci;
            ordering.mapIndex(di.view(),ci);
            values[k] = globalSol[ci];
        }

        // local function space for solid displacement
        typedef typename LocalFunctionSpace::template Child<0>::Type ScalarDispLFS;
        // further types required for gradient calculations
        typedef typename ScalarDispLFS::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
        typedef typename ScalarDispLFS::Traits::FiniteElementType::
                Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef Dune::FieldMatrix<RF, dim, dim> Tensor;

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        (*this)[scvIdx].effPorosity = 0.0;

        auto dispSize = localFunctionSpace.child(0).size();

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
        {
            GlobalPosition scvCenter = fvGeometry.subContVol[scvIdx].localCenter;
            GlobalPosition scvCenterGlobal = element.geometry().global(fvGeometry.subContVol[scvIdx].localCenter);

            // evaluate gradient of displacement shape functions at the center of
            // the sub control volume in the reference element
            std::vector<JacobianType_V> vRefShapeGradient(dispSize);
            localFunctionSpace.child(0).finiteElement().localBasis().evaluateJacobian(scvCenter, vRefShapeGradient);

            // transform gradient to element in global coordinates
            const JacobianInverseTransposed jacInvT = element.geometry().jacobianInverseTransposed(scvCenter);
            std::vector<Dune::FieldVector<RF,dim> > vShapeGradient(dispSize);

            // loop over element vertices
            for (size_t i = 0; i < dispSize; i++)
            {
                    vShapeGradient[i] = 0.0;
                    jacInvT.umv(vRefShapeGradient[i][0],vShapeGradient[i]);
            }

            // calculate gradient of current displacement
            // (gradient of a vector is a tensor)
            Tensor uGradient(0.0);

            // loop over coordinate directions
            for(int coordDir = 0; coordDir < dim; ++coordDir)
            {
                const ScalarDispLFS & scalarDispLFS = localFunctionSpace.child(coordDir);
                // loop over element vertices
                for (size_t i = 0; i < scalarDispLFS.size(); i++)
                    uGradient[coordDir].axpy((*this)[i].displacement(coordDir), vShapeGradient[i]);
            }

            (*this)[scvIdx].volumetricStrain = 0.0;

            for(int coordDir = 0; coordDir < dim; ++coordDir)
            {
                (*this)[scvIdx].volumetricStrain += uGradient[coordDir][coordDir];
            }

            // calculate the effective porosity
            if(problem.coupled() == true)
            {
//                     if ((*this)[scvIdx].divU < -(*this)[scvIdx].initialPorosity())
//                     {
//                         (*this)[scvIdx].effPorosity = (*this)[scvIdx].initialPorosity();
//                         std::cout<<"volume change too large"<<std::endl;
//
//                     }
//                 else
//                     // this equation would be correct if the bulk volume could change (Vol_new = Vol_   *(1+div u)), however, we
//                     // have a constant bulk volume therefore we should apply phi_eff = phi_init + div u
//                     // but this causes convergence problems. Since div u is very small here the chosen relation is
//                     // assumed to be a good approximation
//                      (*this)[scvIdx].effPorosity = ((*this)[scvIdx].initialPorosity() + (*this)[scvIdx].divU)/(1.0 + (*this)[scvIdx].divU);
// //                      std::cout << "effPorosity is " << (*this)[scvIdx].effPorosity << std::endl;
//
                    (*this)[scvIdx].effPorosity = ((*this)[scvIdx].initialPorosity() + (*this)[scvIdx].volumetricStrain)/*/(1.0 + (*this)[scvIdx].volumetricStrain)*/;
                    // Alternative after Cappa & Rutquist (2011)
//                    (*this)[scvIdx].effPorosity = 1 - (1 - (*this)[scvIdx].initialPorosity() )*exp( -((*this)[scvIdx].volumetricStrain) );
//                    std::cout << "effPorosity[" << scvIdx << "] is " << (*this)[scvIdx].effPorosity << std::endl;
            }
            else
                (*this)[scvIdx].effPorosity = (*this)[scvIdx].initialPorosity();

            if( ((scvCenterGlobal[0] > 740) && (scvCenterGlobal[0] < 1260)) &&
            ((scvCenterGlobal[1] > 740) && (scvCenterGlobal[1] < 1260)) )
            {
//                 std::cout << "volumetricStrain[" << eIdx << "][" << scvIdx << "] = " << (*this)[scvIdx].volumetricStrain << std::endl;
//
//                 std::cout << "effPorosity_ = " << (*this)[scvIdx].effPorosity << std::endl;
            }

        }
    }


    const std::vector<Scalar>& dofValues() const
    {
        return dofValues_;
    }

    Scalar& dofValues(int k)
    {
        return dofValues_[k];
    }

    const std::vector<Scalar>& prevValues() const
    {
        return prevValues_;
    }

private:
    std::vector<Scalar> dofValues_;
    std::vector<Scalar> prevValues_;
};

} // namespace Dumux

#endif
