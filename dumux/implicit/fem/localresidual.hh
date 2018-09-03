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
 * \brief Calculates the element-wise residual of models based on the finite element method.
 */
#ifndef DUMUX_FEM_LOCAL_RESIDUAL_HH
#define DUMUX_FEM_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include "properties.hh"

#include <dune/istl/io.hh>

namespace Dumux
{
/*!
 * \ingroup FemModel
 * \ingroup FemLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit finite element method.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class FemLocalResidual
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);

    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FEBasis = typename GET_PROP_TYPE(TypeTag, FeBasis);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int qOrder = GET_PROP_VALUE(TypeTag, FemQuadratureOrder);

    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementGeometry = typename Element::Geometry;
    using ReferenceElements = typename Dune::ReferenceElements<ctype, dim>;

    using LocalView = typename FEBasis::LocalView;
    using LocalIndexSet = typename FEBasis::LocalIndexSet;

public:
    using FluxTermType = Dune::FieldMatrix<Scalar, numEq, dimWorld>;

    //added afterwards
    using DimVector = Dune::FieldVector<Scalar, dim>;


    // copying the local residual class is not a good idea
    FemLocalResidual(const FemLocalResidual&) = delete;

    // the default constructor
    FemLocalResidual() = default;

    /*!
     * \brief Initialize the local residual.
     *
     * This assumes that all objects of the simulation have been fully
     * allocated but not necessarily initialized completely.
     *
     * \param problem The representation of the physical problem to be
     *             solved.
     * \param feBasis The global finite element basis.
     */
    void init(const Problem& problem)
    {
        problemPtr_ = &problem;
        feBasisPtr_ = &problem.model().jacobianAssembler().feBasis();
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    void eval(const Element &element)
    {
        // prepare the element solutions etc...
        const auto& curSol = model().curSol();
        const auto& prevSol = model().prevSol();

        // prepare the current and previous element solutions
        auto localView = feBasis().localView();
        auto localIndexSet = feBasis().localIndexSet();
        localView.bind(element);
        localIndexSet.bind(localView);

        auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
        ElementSolutionVector curElemSol(numLocalDofs);
        ElementSolutionVector prevElemSol(numLocalDofs);
        for (int i = 0; i < numLocalDofs; i++)
        {
            auto dofIdxGlobal = localIndexSet.index(i);
            curElemSol[i] = curSol[dofIdxGlobal];
            prevElemSol[i] = prevSol[dofIdxGlobal];
        }

       // printvector(std::cout, curElemSol, "curElemSol", "");

        // call the eval routine using the prepared local variables
        asImp_().eval(element, localView, localIndexSet, curElemSol, prevElemSol);
    }

    /*!
     * \brief Compute the integral of the storage term in an element.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    PrimaryVariables evalStorage(const Element &element)
    {   std::cout<<"evalStorage wird aufgerufen"<<std::endl;
        // prepare the element solutions etc...
        const auto& curSol = model().curSol();

        // prepare the current and previous element solutions
        auto localView = feBasis().localView();
        auto localIndexSet = feBasis().localIndexSet();
        localView.bind(element);
        localIndexSet.bind(localView);

        auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
        ElementSolutionVector curElemSol(numLocalDofs);
        for (int i = 0; i < numLocalDofs; i++)
            curElemSol[i] = curSol[localIndexSet.index(i)];

        // call the eval routine using the prepared local variables
        return asImp_().evalStorage_(element, localView, localIndexSet, curElemSol);
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual ought to be calculated
     * \param localView The finite element basis bound to the actual element
     * \param localIndexSet The index set bound to the actual element
     * \param curElemSol The current solution at the dofs connected to the element
     * \param prevElemSol The previous solution at the dofs connected to the element
     * \param bcTypes The types of the boundary conditions for all the vertices of the element
     */
    void eval(const Element& element,
              const LocalView& localView,
              const LocalIndexSet& localIndexSet,
              const ElementSolutionVector& curElemSol,
              const ElementSolutionVector& prevElemSol)
    {
        // resize the vectors for all terms
        auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
        residual_.resize(numLocalDofs);
        storageTerm_.resize(numLocalDofs);

        residual_ = 0.0;
        storageTerm_ = 0.0;

       // printvector(std::cout, residual_, "residual", "");

        // integrate over the element and eventually handle boundary conditions
        auto geometry = element.geometry();
        evalVolume_(element, geometry, localView, localIndexSet, curElemSol, prevElemSol);
        evalBoundary_(element, geometry, localView, localIndexSet, curElemSol, prevElemSol);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param element The finite element
     * \param ipData Data on shape values and gradients at the integration point
     * \param secVars Secondary variables and parameters of the problem
     *
     */
    PrimaryVariables computeSource(const Element& element,
                                   const IpData& ipData,
                                   const SecondaryVariables& secVars,
                                   const ElementSolutionVector& elemSol) const
    {
        PrimaryVariables source(0);
//std::cout<<"fem/locres:computeSource wird aufgerufen"<<std::endl;
        // add contributions from volume flux sources
        source += this->problem().source(element, ipData, secVars);

        // TODO: add contribution from possible point sources
        // source += this->problem().scvPointSources(element, fvGeometry, elemVolVars, scv);

        return source;
    }

    /*!
     * \brief Calculate stabilization terms
     *
     * \param element The finite element
     * \param ipData Data on shape values and gradients at the integration point
     * \param secVars Secondary variables and parameters of the problem
     *
     */
    PrimaryVariables computeStabilizationTerms(const Element& element,
                                               const IpData& ipData,
                                               const SecondaryVariables& secVars,
                                               const ElementSolutionVector& elemSol) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Return the problem we are solving. Only call this after init()!
     */
    const Problem& problem() const
    { return *problemPtr_; }

    /*!
     * \brief Return the model used for the simulation. Only call this after init()!
     */
    const Model& model() const
    { return problem().model(); }

    /*!
     * \brief Return the global finite element basis. Only call this after init()!
     */
    const FEBasis& feBasis() const
    { return *feBasisPtr_; }

    /*!
     * \brief Returns the local residual for all dofs of the element.
     */
    const ElementSolutionVector& residual() const
    { return residual_; }

    /*!
     * \brief Returns the local residual for a given local dof inside the element.
     *
     * \param localDofIdx The local dof
     */
    const PrimaryVariables& residual(const int localDofIdx) const
    { return residual_[localDofIdx]; }

protected:

    PrimaryVariables evalStorage_(const Element& element,
                                  const LocalView& localView,
                                  const LocalIndexSet& localIndexSet,
                                  const ElementSolutionVector& curElemSol)
    {
        const auto& localBasis = localView().tree().finiteElement().localBasis();
        auto geometry = element.geometry();

        // initialize container
        PrimaryVariables result(0.0);

        // select quadrature rule
        const auto& rule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), qOrder);

        // loop over quadrature points
        for (auto it = rule.begin(); it != rule.end(); ++it)
        {
            // Obtain and store shape function values and gradients at the current quad point
            IpData ipData(geometry, it->position(), localBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            SecondaryVariables curSecVars;
            curSecVars.update(curElemSol, problem(), element, ipData);

            // evaluate storage term contribution
            PrimaryVariables storage = asImp_().computeStorage(element, ipData, curSecVars);

            // scale by extrusion factors
            storage *= curSecVars.extrusionFactor();

           // std::cout << "curSecVars.extrusionfactor: " << curSecVars.extrusionFactor() << std::endl;


            // Scale by determinant of the transformation and integration weight
            storage *= it->weight() * geometry.integrationElement(it->position());

        //std::cout << "hello" << std::endl;
        //printvector(std::cout, storage, "storage","");

            // add to container
            result += storage;
        }

        return result;
    }

    void evalVolume_(const Element& element,
                     const ElementGeometry& geometry,
                     const LocalView& localView,
                     const LocalIndexSet& localIndexSet,
                     const ElementSolutionVector& curElemSol,
                     const ElementSolutionVector& prevElemSol)
    {
        const auto& localBasis = localView.tree().finiteElement().localBasis();
        auto numLocalDofs = localBasis.size();

        // select quadrature rule
        const auto& rule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), qOrder);
int count = 0;
        // loop over quadrature points
        for (auto it = rule.begin(); it != rule.end(); ++it)
        {   count++;
//std::cout << "femLocResItRuleCounter" << count << std::endl;
            // Obtain and store shape function values and gradients at the current quad point
            IpData ipData(geometry, it->position(), localBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            SecondaryVariables curSecVars, prevSecVars;

//std::cout << "femLocRes vor curSecVars und prevSecVars.update()" << std::endl;
        //printvector(std::cout, curElemSol, "femLocalResCurElemSol","");
        //printvector(std::cout, prevElemSol, "femLocalResPrevElemSol","");

//std::cout << "femLocResCurElemSol: " << std::endl;
            curSecVars.update(curElemSol, problem(), element, ipData);
//std::cout << "femLocResPrevElemSol: " << std::endl;
            prevSecVars.update(prevElemSol, problem(), element, ipData);

//printvector(std::cout, curElemSol, "curElemSol", "");
//    printvector(std::cout, curSecVars.velocity(), "curSecVarsVelocity", "");
//    std::cout << "end" << std::endl;
//printvector(std::cout, prevElemSol, "prevElemSol", "");

            // evaluate storage term contribution
            PrimaryVariables storage = asImp_().computeStorage(element, ipData, curSecVars, curElemSol);
            PrimaryVariables prevStorage = asImp_().computeStorage(element, ipData, prevSecVars, prevElemSol);

            // evaluate source term contribution
            PrimaryVariables source = asImp_().computeSource(element, ipData, curSecVars, curElemSol);


            // evaluate flux term contribution
            FluxTermType flux = asImp_().computeFlux(element, ipData, curSecVars, curElemSol);

            // evaluate stabilization term contributions
//            PrimaryVariables stabTerms = asImp_().computeStabilizationTerms(element, ipData, curSecVars, curElemSol);
            DimVector stabTerms = asImp_().computeStabilizationTerms(element, ipData, curSecVars, curElemSol);



//std::cout <<  "femLocResCurSecVarsPressure: " << curSecVars.pressure() << std::endl;
//std::cout <<  "femLocResPrevSecVarsPressure: " << prevSecVars.pressure() << std::endl;

//printvector(std::cout, storage, "femLocResStorage", "");
//printvector(std::cout, prevStorage, "femLocResPrevStorage", "");
//printvector(std::cout, source, "femLocResSourceresidual", "");
//printmatrix(std::cout, flux, "femLocResFlux", "");

      //      std::cout << "dim = " << dim << std::endl;
      //      std::cout << "dimWorld = " << dimWorld << std::endl;

            //TODO: extrusionfactor zugriff
            // scale terms by extrusion factors
        //            storage *= problem().extrusionFactorAtPos(element.geometry().center());
        //            prevStorage *= problem().extrusionFactorAtPos(element.geometry().center());
        //            source *= problem().extrusionFactorAtPos(element.geometry().center());
        //            flux *= problem().extrusionFactorAtPos(element.geometry().center());

            storage *= curSecVars.extrusionFactor();
            prevStorage *= curSecVars.extrusionFactor();
            source *= curSecVars.extrusionFactor();
            flux *= curSecVars.extrusionFactor();

      //      std::cout << "curSecVars.extrusionfactor: " << curSecVars.extrusionFactor() << std::endl;
      //      std::cout << "prevSecVars.extrusionfactor: " << prevSecVars.extrusionFactor() << std::endl;


        //            printvector(std::cout, storage, "storage", "");
        //            printvector(std::cout, prevStorage, "prevStorage", "");
        //            printvector(std::cout, source, "sourceresidual", "");
        //            printmatrix(std::cout, flux, "flux", "");

            // calculate time derivative
            storage -= prevStorage;

//printvector(std::cout, storage, "femLocResStorageDiff", "");


            storage /= problem().timeManager().timeStepSize();

   //         std::cout << "sind in evalVolume in fem/locres" <<std::endl;
   //		    std::cout << residual_[0][0] << std::endl;

            // add entries to residual vector
            Scalar qWeight = it->weight() * geometry.integrationElement(it->position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                for (unsigned int i = 0; i < numLocalDofs; ++i)
                {
//            std::cout <<"femLocRes localDof i: "<< i << "    Index eqIdx: " << eqIdx <<std::endl;
//            std::cout <<"femLocRes curSecVarsPressure: " << curSecVars.pressure() << std::endl;
//            std::cout <<"femLocRes prevSecVarsPressure: " << prevSecVars.pressure() << std::endl;

        //  std::cout <<"LocResstorage[eqIdx]: "<< storage[eqIdx] <<std::endl;
        //  std::cout <<"LocRessource[eqIdx]: "<< source[eqIdx] <<std::endl;
        //  std::cout <<"LocResflux[eqIdx]: "<< flux[eqIdx] <<std::endl;
        //  printmatrix(std::cout, flux, "flux", "");
        //  std::cout <<"stabTerms[eqIdx]: "<< stabTerms[eqIdx] <<std::endl;
        //  std::cout <<"FemLocResqWeight: "<< qWeight <<std::endl;
                    residual_[i][eqIdx] += (storage[eqIdx] - source[eqIdx])*ipData.shapeValues(i)*qWeight;
                    residual_[i][eqIdx] -= (flux[eqIdx]*ipData.shapeGradients(i))*qWeight;
//		  std::cout <<"FemLocResResidual_[i][eqIdx]: "<< (storage[eqIdx] - source[eqIdx])*ipData.shapeValues(i)*qWeight- (flux[eqIdx]*ipData.shapeGradients(i))*qWeight <<std::endl;
//		  std::cout <<"FemLocResStabTerms: "<<  (stabTerms*ipData.shapeGradients(i))*qWeight<<std::endl;

//                    residual_[i][eqIdx] += stabTerms[eqIdx];
                    if(eqIdx <= 1){
                    residual_[i][eqIdx] += (stabTerms[eqIdx]*curSecVars.velocity()[eqIdx]*ipData.shapeGradients(i)[eqIdx])*qWeight;
                    }
//std::cout <<"FemLocResResidual_[i][eqIdx]: "<< (storage[eqIdx] - source[eqIdx])*ipData.shapeValues(i)*qWeight- (flux[eqIdx]*ipData.shapeGradients(i))*qWeight <<std::endl;
//std::cout <<"FemLocResStabTerms: "<<  (stabTerms[eqIdx]*curSecVars.velocity()[eqIdx]*ipData.shapeGradients(i)[eqIdx])*qWeight<<std::endl;
        //  printvector(std::cout, residual_, "LocResresidual", "");
        //std::cout << "FemLocResipData.shapeValues(i): " << ipData.shapeValues(i) << std::endl;
                }
            }

                //printvector(std::cout, residual_, "FemLocResresidual", "");
        }
    }

    void evalBoundary_(const Element& element,
                       const ElementGeometry& geometry,
                       const LocalView& localView,
                       const LocalIndexSet& localIndexSet,
                       const ElementSolutionVector& curElemSol,
                       const ElementSolutionVector& prevElemSol)
    {
        if (element.hasBoundaryIntersections())
            evalNeumann_(element, geometry, localView, localIndexSet, curElemSol);
    }

    void evalNeumann_(const Element& element,
                      const ElementGeometry& geometry,
                      const LocalView& localView,
                      const LocalIndexSet& localIndexSet,
                      const ElementSolutionVector& curElemSol)
    {
        const auto& localBasis = localView.tree().finiteElement().localBasis();
        auto numLocalDofs = localBasis.size();

        // integrate boundary contribution
        for (const auto& is : intersections(problem().gridView(), element))
        {
            // only handle faces on the boundary
            if (!is.boundary())
                continue;

            // only treat faces with neumann boundary conditions
            auto bcTypes = problem().boundaryTypes(element, is);
            assert(!bcTypes.hasOutflow() && "Outflow BCs are not implemented yet for FEM models");
            if (!bcTypes.hasNeumann())
                continue;

            // select quadrature rule for intersection faces (dim-1)
            auto insideGeom = is.geometryInInside();
            const auto& faceRule = Dune::QuadratureRules<Scalar, dim-1>::rule(insideGeom.type(), qOrder);

            // Treat Neumann boundary conditions
            for (auto it = faceRule.begin(); it != faceRule.end(); ++it)
            {
                // position of quadrature point in local and global coordinates of element
                auto local = insideGeom.global(it->position());

                // evaluate basis functions of all all element vertices for quadrature point
                IpData ipData(geometry, local, localBasis);

                // evaluate secondary variables
                SecondaryVariables secVars;
                secVars.update(curElemSol, problem(), element, ipData);

                // evaluate neumann boundary condition
                auto neumannFlux = problem().neumann(element, is, curElemSol, ipData);

                // get quadrature rule weight for intersection
                Scalar qWeight = it->weight();
                qWeight *= is.geometry().integrationElement(it->position());
                qWeight *= secVars.extrusionFactor();

                // add entries to residual vector
                for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    for (unsigned int i = 0; i < numLocalDofs; ++i)
                        if (bcTypes.isNeumann(eqIdx))
                            residual_[i][eqIdx] += ipData.shapeValues(i)*qWeight*neumannFlux[eqIdx];
            }
        }
//std::cout << "impLocresEvalNeumann " << std::endl;
    }

  private:
    const Implementation& asImp_() const
    { return static_cast<const Implementation&>(*this); }

    Implementation& asImp_()
    { return static_cast<Implementation&>(*this); }

    const Problem* problemPtr_;
    const FEBasis* feBasisPtr_;

    ElementSolutionVector storageTerm_;
    ElementSolutionVector residual_;
};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
