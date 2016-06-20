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
 * \brief This file contains the data which is required to calculate
 *        the mechanic stresses according to Hooke's law.
 */
#ifndef DUMUX_GEOMECHANICS_HOOKES_LAW_HH
#define DUMUX_GEOMECHANICS_HOOKES_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
}

/*!
 * \ingroup CCTpfaHookesLaw
 * \brief Evaluates the stresses, tractions and compressions on a face according to Hooke's law.
 *        Specializations are given for the different discretization methods.
 */
template <class TypeTag, typename DiscretizationMethod = void>
class HookesLaw
{};

template <class TypeTag>
class HookesLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa>::type >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet::IndexType IndexType;
    typedef typename std::vector<IndexType> Stencil;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension} ;
    static constexpr int voigtDim = 0.5*(dim*dim+dim);

    typedef Dune::FieldMatrix<Scalar, voigtDim, voigtDim> StiffnessMatrix;
    typedef Dune::FieldVector<Scalar, voigtDim> VoigtVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    struct FaceData
    {
        Scalar insideLambda, insideMu;
        DimVector insideAlpha, insideN, insideU;

        Scalar outsideLambda, outsideMu;
        DimVector outsideAlpha, outsideN, outsideU;

        bool valueSet;

        FaceData()
        {
            valueSet = false;
        }
    };

public:

    static DimVector stressVector(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        DimMatrix sigma = calculateSigma_(problem, scvFace);

        // calculate Sigma*n
        DimVector stressVec(0.0);
        sigma.mv(scvFace.unitOuterNormal(), stressVec);
        stressVec *= scvFace.area();

        return stressVec;
    }

    static DimMatrix stressTensor(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        return calculateSigma_(problem, scvFace);
    }

    static Stencil stencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        std::vector<IndexType> stencil;
        if (!scvFace.boundary())
        {
            stencil.push_back(scvFace.insideScvIdx());
            stencil.push_back(scvFace.outsideScvIdx());
        }
        else
            stencil.push_back(scvFace.insideScvIdx());

        return stencil;
    }

    static DimMatrix calculateInversA(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        FaceData faceData = obtainFaceData_(problem, scvFace);

        DimMatrix inversA(0.0);
        addEntriesToMatrix_(faceData.insideLambda, faceData.insideMu, faceData.insideAlpha, faceData.insideN, inversA);
        addEntriesToMatrix_(faceData.outsideLambda, faceData.outsideMu, faceData.outsideAlpha, faceData.outsideN, inversA);
        inversA.invert();

        return inversA;
    }

    static DimVector interpolateFaceDisplacement(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        FaceData faceData = obtainFaceData_(problem, scvFace);
        return interpolateFaceDisplacement_(problem, scvFace, faceData);
    }

private:

    static FaceData obtainFaceData_(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        FaceData container;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto& insideVolVars = problem.model().curVolVars(insideScvIdx);
        container.insideU = insideVolVars.displacement();
        container.insideLambda = insideVolVars.lambda();
        container.insideMu = insideVolVars.mu();
        container.insideN = scvFace.unitOuterNormal();
        container.insideAlpha = scvFace.center();
        container.insideAlpha -= insideScv.center();
        container.insideAlpha /= container.insideAlpha.two_norm2();

        const auto outsideScvIdx = scvFace.outsideScvIdx();
        const auto& outsideVolVars = problem.model().curVolVars(outsideScvIdx);
        container.outsideU = outsideVolVars.displacement();
        container.outsideLambda = outsideVolVars.lambda();
        container.outsideMu = outsideVolVars.mu();
        container.outsideN = scvFace.unitOuterNormal();
        container.outsideN *= -1;
        if (scvFace.boundary())
            container.outsideAlpha = 0.0;
        else
        {
            const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
            container.outsideAlpha = scvFace.center();
            container.outsideAlpha -= outsideScv.center();
            container.outsideAlpha /= container.outsideAlpha.two_norm2();
        }

        container.valueSet = true;

        return container;
    }

    static DimMatrix calculateSigma_(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        DimMatrix sigma(0.0);
        StiffnessMatrix C(0.0);
        VoigtVector voigtStrain(0.0);
        VoigtVector voigtSigma(0.0);

        FaceData faceData = obtainFaceData_(problem, scvFace);
        DimVector faceU = interpolateFaceDisplacement_(problem, scvFace, faceData);

        fillStiffnessMatrix_(C, faceData.insideLambda, faceData.insideMu);
        fillStrainVector_(voigtStrain, faceData.insideAlpha, faceData.insideU, faceU);

        C.mv(voigtStrain, voigtSigma);

        if (dim == 2)
        {
            sigma[0][0] = voigtSigma[0];
            sigma[0][1] = voigtSigma[2];
            sigma[1][0] = voigtSigma[2];
            sigma[1][1] = voigtSigma[1];
        }
        else
            DUNE_THROW(Dune::NotImplemented, "dim = " << dim << " is not implemented yet");

        return sigma;
    }

    static DimVector interpolateFaceDisplacement_(const Problem& problem, const SubControlVolumeFace& scvFace, const FaceData& faceData, const bool oldSol = false)
    {
        DimVector faceU(0.0);

        if (!scvFace.boundary())
        {
            DimMatrix inversA(0.0);
            DimMatrix insideB(0.0);
            DimMatrix outsideB(0.0);

            getInversA_(problem, scvFace, faceData, inversA);
            addEntriesToMatrix_(faceData.insideLambda, faceData.insideMu, faceData.insideAlpha, faceData.insideN, insideB);
            addEntriesToMatrix_(faceData.outsideLambda, faceData.outsideMu, faceData.outsideAlpha, faceData.outsideN, outsideB);

            DimVector insideTmp(0.0);
            DimVector outsideTmp(0.0);
            insideB.mv(faceData.insideU, insideTmp);
            outsideB.mv(faceData.outsideU, outsideTmp);

            insideTmp += outsideTmp;

            inversA.mv(insideTmp, faceU);
        }
        else
        {
            if (!oldSol)
            {
                try { return problem.model().curVolVars(scvFace.outsideScvIdx()).displacement(); }
                catch (Dune::Exception& e)
                {
                    DUNE_THROW(Dune::InvalidStateException, "Error ocurred during the displacement interpolation on a boundary scv face. Only call this method on inner scv faces or pure Dirichlet boundaries with the volvars bound to the element");
                }
            }
            else
            {
                // TODO
                DUNE_THROW(Dune::NotImplemented, "Reconstruction of the previous boundary vol vars not yet implemented");
            }
        }

        return faceU;
    }

    template<typename T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache)>::type getInversA_(const Problem& problem,
                                                                                                  const SubControlVolumeFace& scvFace,
                                                                                                  const FaceData& faceData,
                                                                                                  DimMatrix& inversA)
    { inversA = problem.model().fluxVarsCache(scvFace).inversA(); }

    template<typename T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableFluxVariablesCache)>::type getInversA_(const Problem& problem,
                                                                                                   const SubControlVolumeFace& scvFace,
                                                                                                   const FaceData& faceData,
                                                                                                   DimMatrix& inversA)
    {
        addEntriesToMatrix_(faceData.insideLambda, faceData.insideMu, faceData.insideAlpha, faceData.insideN, inversA);
        addEntriesToMatrix_(faceData.outsideLambda, faceData.outsideMu, faceData.outsideAlpha, faceData.outsideN, inversA);
        inversA.invert();
    }

    static void addEntriesToMatrix_(const Scalar lambda, const Scalar mu, const DimVector& alpha, const DimVector& normal, DimMatrix& matrix)
    {
        if (dim == 2)
        {
            matrix[0][0] += (lambda + 2*mu)*alpha[0]*normal[0] + mu*alpha[1]*normal[1];
            matrix[0][1] += lambda*alpha[1]*normal[0] + mu*alpha[0]*normal[1];
            matrix[1][0] += mu*alpha[1]*normal[0] + lambda*alpha[0]*normal[1];
            matrix[1][1] += mu*alpha[0]*normal[0] + (lambda + 2*mu)*alpha[1]*normal[1];
        }
        else
            DUNE_THROW(Dune::NotImplemented, "dim = " << dim << " is not implemented yet");
    }

    static void fillStiffnessMatrix_(StiffnessMatrix& C, const Scalar lambda, const Scalar mu)
    {
        if (dim == 2)
        {
            C[0][0] = lambda + 2*mu;
            C[0][1] = lambda;
            C[0][2] = 0.0;

            C[1][0] = lambda;
            C[1][1] = lambda + 2*mu;
            C[1][2] = 0.0;

            C[2][0] = 0.0;
            C[2][1] = 0.0;
            C[2][2] = mu;
        }
    }

    static void fillStrainVector_(VoigtVector& strain, const DimVector& alpha, const DimVector& insideU, const DimVector& faceU)
    {
        if (dim == 2)
        {
            strain[0] = alpha[0]*(faceU[0] - insideU[0]);
            strain[1] = alpha[1]*(faceU[1] - insideU[1]);
            strain[2] = alpha[1]*(faceU[0] - insideU[0]) + alpha[0]*(faceU[1] - insideU[1]);
        }
    }
};

} // end namespace

#endif
