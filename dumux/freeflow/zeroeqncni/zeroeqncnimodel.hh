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
 * \brief Base class for all models which use the non-isothermal
 *        compositional ZeroEq box model.
 */
#ifndef DUMUX_ZEROEQNCNI_MODEL_HH
#define DUMUX_ZEROEQNCNI_MODEL_HH

#include "zeroeqncniindices.hh"
#include "zeroeqncniproperties.hh"
#include <dumux/freeflow/zeroeqnc/zeroeqncmodel.hh>


namespace Dumux
{
/*!
 * \ingroup BoxZeroEqncniModel
 * \brief Adaptation of the box scheme to the non-isothermal compositional ZeroEq model.
 *
 * This model implements an single-phase non-isothermal compositional free flow
 * solving the mass and the momentum balance. For the momentum balance
 * the Reynolds-averaged Navier-Stokes (RANS) equation with zero equation
 * (algebraic) turbulence model is used.
 *
 * Mass balance:
 * \f[
 *  \frac{\partial \varrho_\textrm{g}}{\partial t}
 *  + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g}\right)
 *  - q_\textrm{g} = 0
 * \f]
 *
 * Momentum Balance:
 * \f[
 *   \frac{\partial \left(\varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g}\right)}{\partial t}
 *   + \boldsymbol{\nabla} \boldsymbol{\cdot} \left(
 *     \varrho_\textrm{g} {\boldsymbol{v}_\textrm{g} {\boldsymbol{v}}_\textrm{g}}
 *     - \left[ \mu_\textrm{g} + \mu_\textrm{g,t} \right]
 *       \left(\boldsymbol{\nabla} \boldsymbol{v}_\textrm{g}
 *             + \boldsymbol{\nabla} \boldsymbol{v}_\textrm{g}^T \right)
 *   \right)
 *   + \left(p_\textrm{g} {\bf {I}} \right)
 *   - \varrho_\textrm{g} {\bf g} = 0,
 * \f]
 *
 * Component mass balance equations:
 * \f[
 *  \frac{\partial \left(\varrho_\textrm{g} X_\textrm{g}^\kappa\right)}{\partial t}
 *  + \boldsymbol{\nabla} \boldsymbol{\cdot} \left( \varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g} X_\textrm{g}^\kappa
 *  - \left[ D^\kappa_\textrm{g} + D^\kappa_\textrm{g,t} \right]
 *    \varrho_\textrm{g} \frac{M^\kappa}{M_\textrm{g}} \boldsymbol{\nabla} x_\textrm{g}^\kappa \right)
 *  - q_\textrm{g}^\kappa = 0
 * \f]
 *
 * Energy balance equation:
 * \f[
 *  \frac{\partial (\varrho_\textrm{g}  u_\textrm{g})}{\partial t}
 *  + \boldsymbol{\nabla} \boldsymbol{\cdot} \left( \varrho_\textrm{g} h_\textrm{g} {\boldsymbol{v}}_\textrm{g}
 *  - \sum_\kappa \left( h^\kappa_\textrm{g} \left[ D^\kappa_\textrm{g} + D^\kappa_\textrm{g,t} \right]
 *                       \varrho_\textrm{g} \frac{M^\kappa}{M_\textrm{g}} \nabla x^\kappa_\textrm{g} \right)
 *  - \left[ \lambda_\textrm{g} + \lambda_\textrm{g,t} \right] \boldsymbol{\nabla} T \right)
 *  - q_\textrm{T} = 0
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class ZeroEqncniModel : public ZeroEqncModel<TypeTag>
{
    typedef ZeroEqncModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        intervals = GET_PROP_VALUE(TypeTag, NumberOfIntervals),
        walls = (GET_PROP_VALUE(TypeTag, BBoxMinIsWall) ? 1 : 0)
                + (GET_PROP_VALUE(TypeTag, BBoxMaxIsWall) ? 1 : 0),
        prec = Indices::scvDataPrecision, // precision of scv data
        width = Indices::scvDataWidth // width of column
    };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { transportCompIdx = Indices::transportCompIdx,
           numComponents = Indices::numComponents };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;


public:
    ZeroEqncniModel()
        : flowNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, FlowNormal))
        , wallNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, WallNormal))
    {
        eps_ = 1e-6;

        // check whether sand grain roughness may be used
        if (GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMinSandGrainRoughness) > 0
            || GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMaxSandGrainRoughness) > 0)
        {
            Dune::dwarn << "warning: surface roughness will not be used for eddy conductivity models." << std::endl;
        }
    }

    //! \copydoc ImplicitModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        const Scalar scale_ = GET_PROP_VALUE(TypeTag, Scaling);

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);

        ScalarField *moleFraction[numComponents];
        for (int i = 0; i < numComponents; ++i)
            moleFraction[i] = writer.template allocateManagedBuffer<Scalar, 1>(numVertices);
        ScalarField *massFraction[numComponents];
        for (int i = 0; i < numComponents; ++i)
            massFraction[i] = writer.template allocateManagedBuffer<Scalar, 1>(numVertices);

        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);
        ScalarField &T = *writer.allocateManagedBuffer(numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            int idx = this->elementMapper().map(*eIt);
            rank[idx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            int numLocalVerts = eIt->template count<dim>();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int vIdxGlobal = this->vertexMapper().map(*eIt, i, dim);
                volVars.update(sol[vIdxGlobal],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               i,
                               false);

                pN[vIdxGlobal] = volVars.pressure()*scale_;
                delP[vIdxGlobal] = volVars.pressure()*scale_ - 1e5;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    (*moleFraction[compIdx])[vIdxGlobal]= volVars.moleFraction(compIdx);
                    (*massFraction[compIdx])[vIdxGlobal]= volVars.massFraction(compIdx);
                    Valgrind::CheckDefined((*moleFraction[compIdx])[vIdxGlobal]);
                    Valgrind::CheckDefined((*massFraction[compIdx])[vIdxGlobal]);
                }
                rho[vIdxGlobal] = volVars.density()*scale_*scale_*scale_;
                mu[vIdxGlobal] = volVars.dynamicViscosity()*scale_;
                velocity[vIdxGlobal] = volVars.velocity();
                velocity[vIdxGlobal] *= 1/scale_;
                T[vIdxGlobal] = volVars.temperature();
            }
        }
        writer.attachVertexData(pN, "P");
        writer.attachVertexData(delP, "delP");

        for (int j = 0; j < numComponents; ++j)
        {
            std::ostringstream moleFrac, massFrac;
            moleFrac << "x_" << FluidSystem::phaseName(phaseIdx)
                     << "^" << FluidSystem::componentName(j);
            writer.attachVertexData(*moleFraction[j], moleFrac.str().c_str());

            massFrac << "X_" << FluidSystem::phaseName(phaseIdx)
                     << "^" << FluidSystem::componentName(j);
            writer.attachVertexData(*massFraction[j], massFrac.str().c_str());
        }

        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);
        writer.attachVertexData(T, "temperature");


        // just to have the actual values plotted
        asImp_().updateWallProperties();

        eIt = this->gridView_().template begin<0>();
        eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false);

            IntersectionIterator isIt = this->gridView_().ibegin(*eIt);
            const IntersectionIterator &endIt = this->gridView_().iend(*eIt);

            for (; isIt != endIt; ++isIt)
            {
                int fIdx = isIt->indexInInside();

                FluxVariables fluxVars(this->problem_(),
                                                    *eIt,
                                                    fvGeometry,
                                                    fIdx,
                                                    elemVolVars,
                                                    false);

                GlobalPosition globalPos = fvGeometry.subContVolFace[fIdx].ipGlobal;

                if (asImp_().shouldWriteSCVData(globalPos))
                asImp_().writeSCVData(volVars, fluxVars, globalPos);
            }
        }
    }

    //! \copydoc ZeroEqModel::writeSCVHeader
    void writeSCVHeader(std::stringstream &stream, const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        ParentType::writeSCVHeader(stream, fluxVars, globalPos);

        stream << std::setw(width) << "eddyCMod"
               << std::setw(width) << GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyConductivityModel)
               << " - " << std::setw(width*2-3) << std::left  << eddyConductivityModelName() << std::right;
    }

    //! \copydoc ZeroEqModel::writeDataHeader
    void writeDataHeader(std::stringstream &stream, int posIdx)
    {
        ParentType::writeDataHeader(stream, posIdx);

        stream << std::setw(width) << "T [K]"
               << std::setw(width) << "T/T_max [-]"
               << std::setw(width) << "TGrad [K/m]"
               << std::setw(width) << "lambda [W/(m*K)]"
               << std::setw(width) << "lambda_t [W/(m*K)]"
               << std::setw(width) << "Pr [-]"
               << std::setw(width) << "Pr_t [-]";
    }

    //! \copydoc ZeroEqModel::writeSCVDataValues
    void writeSCVDataValues(std::stringstream &stream, const VolumeVariables &volVars, const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        ParentType::writeSCVDataValues(stream, volVars, fluxVars, globalPos);

        int posIdx = this->getPosIdx(globalPos);
        int wallIdx = this->getWallIdx(globalPos, posIdx);
        stream << std::setprecision(prec)
               << std::setw(width) << fluxVars.temperature()
               << std::setw(width) << fluxVars.temperature() / this->wall[wallIdx].maxTemperature[posIdx]
               << std::setw(width) << fluxVars.temperatureGrad()[wallNormal_]
               << std::setw(width) << fluxVars.thermalConductivity()
               << std::setw(width) << fluxVars.thermalEddyConductivity()
               << std::setw(width) << fluxVars.prandtlNumber()
               << std::setw(width) << fluxVars.turbulentPrandtlNumber();
    }


    //! \copydoc ZeroEqModel::resetWallProperties
    void resetWallProperties()
    {
        ParentType::resetWallProperties();

        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
            for (int posIdx = 0; posIdx < intervals; ++posIdx)
            {
                this->wall[wallIdx].maxTemperature[posIdx] = 0.0;
            }
    }

    //! \copydoc ZeroEqModel::calculateMaxFluxVars
    void calculateMaxFluxVars(const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        ParentType::calculateMaxFluxVars(fluxVars, globalPos);

        int posIdx = this->getPosIdx(globalPos);
        int wallIdx = this->getWallIdx(globalPos, posIdx);
        if (this->wall[wallIdx].maxTemperature[posIdx] < fluxVars.temperature())
            for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
                this->wall[wallIdx].maxTemperature[posIdx] = fluxVars.temperature();
    }

    //! \copydoc ZeroEqModel::doInterpolationFluxValues
    const void doInterpolationFluxValues(const int wallIdx, const int posIdx, const int prevIdx, const int nextIdx)
    {
        ParentType::doInterpolationFluxValues(wallIdx, posIdx, prevIdx, nextIdx);
        this->wall[wallIdx].maxTemperature[posIdx] = this->interpolation(posIdx, prevIdx, this->wall[wallIdx].maxTemperature[prevIdx], nextIdx, this->wall[wallIdx].maxTemperature[nextIdx]);
    }

    /*!
     * \brief Returns the name of used the eddy conductivity model.
     */
    const char *eddyConductivityModelName() const
    {
        switch (GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyConductivityModel))
        {
            case EddyConductivityIndices::noEddyConductivityModel: // 0
                return "noEddyConductivityModel";
                break;
            case EddyConductivityIndices::reynoldsAnalogy: // 1
                return "reynoldsAnalogy";
                break;
            case EddyConductivityIndices::modifiedVanDriest: // 2
                return "modifiedVanDriest";
                break;
            case EddyConductivityIndices::deissler: // 3
                return "deissler";
                break;
            case EddyConductivityIndices::meier: // 4
                return "meier";
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "This eddy conductivity model is not implemented.");
        }
    }

private:
    Scalar eps_;
    const int flowNormal_;
    const int wallNormal_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#include "zeroeqncnipropertydefaults.hh"

#endif // DUMUX_ZEROEQNCNI_MODEL_HH
