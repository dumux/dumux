// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Orthotropic wood model with moisture-dependent multiplicative growth/shrinkage.
 *
 * Primary variables: displacement \f$ d_i,\; i=0,\ldots,\mathrm{dim}-1 \f$ and
 * moisture content \f$ m \f$.
 *
 * The system of equations is:
 * - Momentum balance (per spatial direction \f$ i \f$):
 *   \f[ -\nabla_X \cdot P = 0 \f]
 *   where \f$ P \f$ is the first Piola-Kirchhoff stress tensor computed from the
 *   multiplicative decomposition \f$ F = F_e F_g \f$.
 *
 * - Moisture transport (pulled back to reference configuration):
 *   \f[ \frac{\partial(J\,m)}{\partial t} - \nabla_X \cdot (J\,F^{-1}\,D\,F^{-T}\,\nabla_X m) = 0 \f]
 */
#ifndef DUMUX_WOOD_TEST_MODEL_HH
#define DUMUX_WOOD_TEST_MODEL_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \brief Equation and primary-variable indices for the wood model.
 * \tparam dim Spatial dimension.
 *
 * Displacement components occupy indices 0 to dim-1;
 * moisture content occupies index dim.
 */
template<int dim>
struct WoodIndices
{
    static constexpr int displacementIdx(int i) { return i; } //!< Index of the i-th displacement primary variable
    static constexpr int momentumEqIdx(int i) { return i; } //!< Index of the i-th momentum balance equation
    static constexpr int moistureIdx = dim; //!< Index of the moisture primary variable
    static constexpr int moistureEqIdx = dim; //!< Index of the moisture transport equation
};

/*!
 * \brief Model traits for the wood model.
 * \tparam dim Spatial dimension.
 */
template<int dim>
struct WoodModelTraits
{
    using Indices = WoodIndices<dim>;
    static constexpr int numEq() { return dim + 1; }
};

/*!
 * \brief Traits class for WoodVolumeVariables.
 * \tparam PV Type of the primary variables vector.
 * \tparam MT The model traits type.
 */
template<class PV, class MT>
struct WoodVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

/*!
 * \brief Volume variables for the wood model.
 *
 * Stores the displacement field \f$ d_i \f$ and the moisture content \f$ m \f$
 * as primary variables.
 *
 * \tparam Traits Traits class providing PrimaryVariables and ModelTraits.
 */
template<class Traits>
class WoodVolumeVariables : public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using Indices = typename Traits::ModelTraits::Indices;

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;

    //! Return the i-th displacement component \f$ d_i \f$
    Scalar displacement(int i) const
    { return this->priVar(Indices::displacementIdx(i)); }

    //! Return the moisture content \f$ m \f$
    Scalar moistureContent() const
    { return this->priVar(Indices::moistureIdx); }
};

/*!
 * \brief Local residual for the wood model.
 *
 * Implements:
 * - Quasi-static momentum balance via the first Piola-Kirchhoff stress tensor
 *   from a multiplicative decomposition \f$ F = F_e F_g \f$ with orthotropic
 *   growth \f$ F_g = Q\,\mathrm{diag}(1+\alpha_R\Delta m,\,1+\alpha_T\Delta m)\,Q^T \f$
 *   and St. Venant-Kirchhoff elasticity in the local material frame.
 * - Pulled-back moisture diffusion:
 *   \f[ \frac{\partial(J\,m)}{\partial t}
 *       - \nabla_X \cdot (J\,F^{-1}\,D\,F^{-T}\,\nabla_X m) = 0 \f]
 *
 * \tparam TypeTag The type tag.
 */
template<class TypeTag>
class WoodLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;
    using ParentType::evalStorage;

    /*!
     * \brief Evaluate storage term (moisture equation only).
     *
     * The pulled-back moisture storage reads:
     * \f[ r_m = \frac{J^{n+1}\,m^{n+1} - J^n\,m^n}{\Delta t}\,V_{\mathrm{scv}} \f]
     * Momentum equations are quasi-static and contribute no storage term.
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        const auto& geom = element.geometry();
        const auto ipLocal = geom.local(scv.dofPosition());
        const auto Fprev = computeF_(geom, ipLocal, fvGeometry, prevElemVolVars);
        const auto Fcur  = computeF_(geom, ipLocal, fvGeometry, curElemVolVars);

        const Scalar mPrev = prevElemVolVars[scv].moistureContent();
        const Scalar mCur  = curElemVolVars[scv].moistureContent();

        const Scalar volume = curElemVolVars[scv].extrusionFactor()
                              * Extrusion::volume(fvGeometry, scv);
        const Scalar dt = this->timeLoop().timeStepSize();

        residual[scv.localDofIndex()][Indices::moistureEqIdx]
            += (Fcur.determinant() * mCur - Fprev.determinant() * mPrev) * volume / dt;
    }

    /*!
     * \brief Evaluate flux terms on a sub-control-volume face.
     *
     * Momentum flux: \f$ -P\,\hat{n}\,\mathrm{d}A \f$, where \f$ P \f$ is the
     * first Piola-Kirchhoff stress tensor.
     *
     * Moisture flux (pulled-back Fick's law on the reference configuration):
     * \f[ q_m = -J\,F^{-1}\,D\,F^{-T}\,\nabla_X m \cdot \hat{N}\,\mathrm{d}A \f]
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        const auto normal = scvf.unitOuterNormal();
        const Scalar area = scvf.area();

        // F = grad(d) + I, and grad(m), m at the integration point
        Tensor F(0.0);
        GlobalPosition gradM(0.0);
        Scalar m = 0.0;
        const auto& shapeValues = fluxVarCache.shapeValues();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            const auto& gradN = fluxVarCache.gradN(scv.indexInElement());
            const Scalar mScv = volVars.moistureContent();
            for (int dir = 0; dir < dimWorld; ++dir)
                F[dir].axpy(volVars.displacement(dir), gradN);
            gradM.axpy(mScv, gradN);
            m += mScv * shapeValues[scv.indexInElement()][0];
        }
        for (int dir = 0; dir < dimWorld; ++dir)
            F[dir][dir] += 1.0;

        NumEqVector flux(0.0);

        // momentum: -P n dA
        const auto P = firstPiolaKirchhoffStressTensor_(problem, element, F, m);
        GlobalPosition Pn(0.0);
        P.mv(normal, Pn);
        for (int dir = 0; dir < dimWorld; ++dir)
            flux[Indices::momentumEqIdx(dir)] = -Pn[dir] * area;

        // moisture: pulled-back Fick on the reference configuration
        //   D (global frame) = Q · diag(D_R, D_T) · Q^T  (already assembled in spatial params)
        //   D_ref = J · F^{-1} · D · F^{-T}              (pulled back to reference config)
        //   flux residual on the SCV face = -vtmv(N, D_ref, ∇_X m) · dA
        const Scalar J = F.determinant();
        Tensor Finv = F; Finv.invert();
        const Tensor FinvT = transpose(Finv);

        Tensor Dref = problem.spatialParams().moistureDiffusivity(scvf.ipGlobal());
        Dref.rightmultiply(FinvT);
        Dref.leftmultiply(Finv);
        Dref *= J;

        flux[Indices::moistureEqIdx]
            = -vtmv(normal, Dref, gradM) * area;

        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        return ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);
    }

private:
    /*!
     * \brief Compute the first Piola-Kirchhoff stress tensor.
     *
     * Uses multiplicative decomposition \f$ F = F_e F_g \f$ with:
     * - Growth tensor \f$ F_g = Q\,\mathrm{diag}(g_R,g_T)\,Q^T \f$,
     *   \f$ g_R = 1+\alpha_R\Delta m \f$, \f$ g_T = 1+\alpha_T\Delta m \f$,
     *   \f$ Q(x) \f$ the spatially-varying rotation from local to global frame.
     * - Elastic deformation gradient \f$ F_e = F\,F_g^{-1} \f$.
     * - Elastic Green-Lagrange strain \f$ E_e = \tfrac{1}{2}(F_e^T F_e - I) \f$.
     * - Orthotropic St. Venant-Kirchhoff law in the local material frame:
     *   \f[ S_e^{\mathrm{loc}} = \begin{pmatrix} C_{11} & C_{12} \\ C_{12} & C_{22} \end{pmatrix}
     *       E_e^{\mathrm{loc}} + 2\,C_{66}\,E_{e,12}^{\mathrm{loc}} \f]
     * - Pull-back to reference: \f$ S = J_g\,F_g^{-1}\,S_e\,F_g^{-T} \f$.
     * - \f$ P = F\,S \f$.
     *
     * Note: \f$ F_g \f$ is symmetric (\f$ F_g^{-T} = F_g^{-1} \f$) but both factors
     * are kept explicit for clarity.
     */
    static Tensor firstPiolaKirchhoffStressTensor_(const Problem& problem,
                                                   const Element& element,
                                                   const Tensor& F,
                                                   Scalar m)
    {
        const auto& sp = problem.spatialParams();
        const Scalar dm = m - sp.referenceMoisture();

        // local rotation Q(x) at the element center (varies smoothly across the blank)
        const auto Q = sp.rotation(element.geometry().center());

        // F_g = Q · diag(1+alpha_R * dm, 1+alpha_T * dm) · Q^T
        const Scalar gR = 1.0 + sp.shrinkageR() * dm;
        const Scalar gT = 1.0 + sp.shrinkageT() * dm;
        const Scalar detFg = gR * gT;

        Tensor FgLocal(0.0);
        FgLocal[0][0] = gR;
        FgLocal[1][1] = gT;

        Tensor FgInvLocal(0.0);
        FgInvLocal[0][0] = 1.0/gR;
        FgInvLocal[1][1] = 1.0/gT;
        const Tensor FgInv = rotateInPlace_(Q, FgInvLocal);

        // F_e = F · F_g^{-1}
        Tensor Fe = F;
        Fe.rightmultiply(FgInv);

        // E_e = 0.5*(F_e^T F_e - I)
        Tensor Ee(0.0);
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                for (int k = 0; k < dimWorld; ++k)
                    Ee[i][j] += Fe[k][i] * Fe[k][j];
        Ee *= 0.5;
        for (int i = 0; i < dimWorld; ++i)
            Ee[i][i] -= 0.5;

        // E_e_local = Q^T · E_e · Q
        const Tensor EeLocal = rotateInverse_(Q, Ee);

        // orthotropic Hooke (plane-stress St. Venant-Kirchhoff) in local frame
        Tensor SeLocal(0.0);
        SeLocal[0][0] = sp.C11() * EeLocal[0][0] + sp.C12() * EeLocal[1][1];
        SeLocal[1][1] = sp.C12() * EeLocal[0][0] + sp.C22() * EeLocal[1][1];
        SeLocal[0][1] = 2.0 * sp.C66() * EeLocal[0][1];
        SeLocal[1][0] = SeLocal[0][1];

        // S_e = Q · S_e_local · Q^T
        const Tensor Se = rotateInPlace_(Q, SeLocal);

        // S = J_g · F_g^{-1} · S_e · F_g^{-T}
        const Tensor FgInvT = transpose(FgInv);

        Tensor S = Se;
        S.rightmultiply(FgInvT);
        S.leftmultiply(FgInv);
        S *= detFg;

        // P = F · S
        Tensor P = F;
        P.rightmultiply(S);
        return P;
    }

    //! Compute \f$ Q\,A\,Q^T \f$: push a tensor from the local material frame to the global frame.
    template<class Rot>
    static Tensor rotateInPlace_(const Rot& Q, const Tensor& A)
    {
        Tensor QA(0.0);
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                for (int k = 0; k < dimWorld; ++k)
                    QA[i][j] += Q[i][k] * A[k][j];
        Tensor R(0.0);
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                for (int k = 0; k < dimWorld; ++k)
                    R[i][j] += QA[i][k] * Q[j][k];
        return R;
    }

    //! Compute \f$ Q^T\,A\,Q \f$: pull a tensor from the global frame to the local material frame.
    template<class Rot>
    static Tensor rotateInverse_(const Rot& Q, const Tensor& A)
    {
        Tensor QtA(0.0);
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                for (int k = 0; k < dimWorld; ++k)
                    QtA[i][j] += Q[k][i] * A[k][j];
        Tensor R(0.0);
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                for (int k = 0; k < dimWorld; ++k)
                    R[i][j] += QtA[i][k] * Q[k][j];
        return R;
    }

    //! Compute the deformation gradient \f$ F = \nabla_X d + I \f$ at a given
    //! local coordinate inside the element from the current element volume variables.
    template<class Geom, class LocalPos, class ElemVV>
    static Tensor computeF_(const Geom& geom,
                            const LocalPos& ipLocal,
                            const FVElementGeometry& fvGeometry,
                            const ElemVV& elemVolVars)
    {
        const auto& localBasis = fvGeometry.feLocalBasis();
        using ShapeJacobian = typename std::decay_t<decltype(localBasis)>::Traits::JacobianType;
        std::vector<ShapeJacobian> shapeJacobian;
        localBasis.evaluateJacobian(ipLocal, shapeJacobian);
        const auto jacInvT = geom.jacobianInverseTransposed(ipLocal);

        std::vector<GlobalPosition> gradN(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            jacInvT.mv(shapeJacobian[scv.localDofIndex()][0], gradN[scv.indexInElement()]);

        Tensor F(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            for (int dir = 0; dir < dimWorld; ++dir)
                F[dir].axpy(volVars.displacement(dir), gradN[scv.indexInElement()]);
        }
        for (int dir = 0; dir < dimWorld; ++dir)
            F[dir][dir] += 1.0;
        return F;
    }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {
struct Wood { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Wood>
{ using type = WoodModelTraits<GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Wood>
{ using type = WoodLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Wood>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = WoodVolumeVariablesTraits<PV, MT>;
    using type = WoodVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
