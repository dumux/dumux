// $Id$
/*****************************************************************************
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_BUCKLEYLEVERETT_ANALYTICAL_HH
#define DUMUX_BUCKLEYLEVERETT_ANALYTICAL_HH

#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

/**
 * @file
 * @brief  Analytical solution of the buckley-leverett problem
 * @author Markus Wolff
 */

namespace Dumux
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * the Buckley-Leverett problem
 */

template<class Scalar, class Law> struct CheckMaterialLaw
{
    static bool isLinear()
    {
        return false;
    }
};

template<class Scalar> struct CheckMaterialLaw<Scalar, LinearMaterial<Scalar> >
{
    static bool isLinear()
    {
        return true;
    }
};

template<class TypeTag> class BuckleyLeverettAnalytic
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum
    {
        dim = GridView::dimension, dimworld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > BlockVector;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:

    // functions needed for analytical solution

    void initializeAnalytic()
    {
        analyticSolution_.resize(size_);
        analyticSolution_ = 0;
        error_.resize(size_);
        error_ = 0;
        elementVolume_.resize(size_);
        elementVolume_ = 0;

        return;
    }

    void prepareAnalytic()
    {
        const MaterialLawParams& materialLawParams(problem_.spatialParameters().materialLawParams(dummyGlobal_, dummyElement_));

        swr_ = materialLawParams.Swr();
        snr_ = materialLawParams.Snr();
        porosity_ = problem_.spatialParameters().porosity(dummyGlobal_, dummyElement_);

        time_ = 0;

        satVec_ = swr_;
        for (int i = 1; i < pointNum_; i++)
        {
            satVec_[i] = satVec_[i - 1] + (1 - snr_ - swr_) / intervalNum_;
        }
//         std::cout<<"satVec_ = "<<satVec_<<std::endl;

        FluidState fluidState;
        Scalar temp = problem_.temperature(dummyGlobal_, dummyElement_);
        Scalar press = problem_.referencePressure(dummyGlobal_, dummyElement_);
        Scalar viscosityW = FluidSystem::phaseViscosity(wPhaseIdx, temp, press, fluidState);
        Scalar viscosityNW = FluidSystem::phaseViscosity(nPhaseIdx, temp, press, fluidState);

        for (int i = 0; i < pointNum_; i++)
        {
            fractionalW_[i] = MaterialLaw::krw(materialLawParams, satVec_[i])/viscosityW;
            fractionalW_[i] /= (fractionalW_[i] + MaterialLaw::krn(materialLawParams, satVec_[i])/viscosityNW);
        }
//         std::cout<<"fractionalW_ = "<<fractionalW_<<std::endl;

        dfwdsw_ = 0;
        for (int i = 1; i < intervalNum_; i++)
        {
            dfwdsw_[i] = (fractionalW_[i + 1] - fractionalW_[i - 1]) / (satVec_[i + 1] - satVec_[i - 1]);
        }
//         std::cout<<"dfwdsw_ = "<<dfwdsw_<<std::endl;

        for (int i = 0; i < pointNum_; i++)
        {
            if (dfwdsw_[i] > dfwdsw_[i + 1])
            {
                dfwdswmax_ = i;
                break;
            }
        }

        return;
    }

    void calcSatError(BlockVector &Approx)
    {
        error_ = 0;
        elementVolume_ = 0;

        ElementIterator eItEnd = problem_.gridView().template end<0> ();
        for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
        {
            // get entity
            const Element& element = *eIt;
            int index = problem_.variables().index(*eIt);
            elementVolume_[index] = element.geometry().volume();
            // std::cout<<"elementVolume_ = "<<elementVolume_[index]<<std::endl;
        }

        Scalar globalVolume = elementVolume_.one_norm();
        // std::cout<<"globalVolume = "<<globalVolume<<std::endl;

        for (int i = 0; i < size_; i++)
        {
            error_[i] = analyticSolution_[i] - Approx[i];
            // std::cout<<"error_ = "<<error_[i]<<std::endl;
            // std::cout<<"analyticSolution_ = "<<analyticSolution_[i]<<std::endl;
            // std::cout<<"Approx = "<<Approx[i]<<std::endl;
        }
        // std::cout<<"error_ = "<<error_<<std::endl;

        Scalar diffNorm = error_.two_norm();
        // std::cout<<"diffNorm = "<<diffNorm<<std::endl;

        for (int i = 0; i < size_; i++)
        {
            error_[i] = diffNorm * pow((elementVolume_[i] / globalVolume), 0.5);
        }

        return;
    }

    void updateExSol()
    {
        // position of the fluid front
        xf_ = 0;
        for (int i = 0; i < pointNum_; i++)
        {
            xf_[i] = vTot_ * time_ / porosity_ * dfwdsw_[i];
        }

//         std::cout<<"xf_ = "<<xf_<<std::endl;
        int xhelp = pointNum_ / 3;
        int xhelpold = 0;
        int xhelpoldold = 0;
        int xfmax = 0;

        // position of maximum xf_
        for (int i = 0; i < pointNum_; i++)
        {
            if (xf_[i] > xf_[i + 1])
            {
                xfmax = i;
                break;
            }
        }

        // balancing of the areas ahead of the front and below the curve
        bool a = true;
        Scalar A1;
        Scalar A2;
        Scalar b;
        int xhelp2 = 0;

        while (a)
        {
            if (CheckMaterialLaw<Scalar, MaterialLaw>::isLinear())
                break;

            A1 = 0;

            for (int i = 0; i <= xhelp - 1; i++)
            {
                A1 += (satVec_[i] - swr_ + satVec_[i + 1] - swr_) * 0.5 * (xf_[i + 1] - xf_[i]);
            }

            A2 = 0;

            for (int i = xhelp; i <= xfmax - 1; i++)
            {
                A2 += (satVec_[xfmax] - satVec_[i] + satVec_[xfmax] - satVec_[i + 1]) * 0.5 * (xf_[i + 1] - xf_[i]);
            }

            b = xf_[xfmax];
            xhelp2 = xfmax;

            while (b > xf_[xhelp])
            {
                xhelp2 += 1;
                b = xf_[xhelp2];
            }

            for (int i = xfmax; i <= xhelp2; i++)
            {
                A2 += (satVec_[i] - satVec_[xfmax] + satVec_[i + 1] - satVec_[xfmax]) * 0.5 * (xf_[i] - xf_[i + 1]);
            }

            xhelpoldold = xhelpold;
            xhelpold = xhelp;

            if (fabs(A1) > fabs(A2))
            {
                xhelp = xhelp - 1;
            }

            if (fabs(A1) < fabs(A2))
            {
                xhelp = xhelp + 1;
            }

            if (xhelp == xhelpoldold)
            {
                a = false;
            }
        }
        // std::cout<<"xf_ = "<<xf_<<std::endl;

        // iterate over vertices and get analytic saturation solution
        ElementIterator eItEnd = problem_.gridView().template end<0> ();
        for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
        {
            // get global coordinate of cell center
            const GlobalPosition& globalPos = eIt->geometry().center();

            // find index of current vertex
            int index = problem_.variables().index(*eIt);

            // account for linear material law
            if (CheckMaterialLaw<Scalar, MaterialLaw>::isLinear())
            {
                if (globalPos[0] <= xf_[1])
                {
                    analyticSolution_[index] = 1 - snr_;
                }
                if (globalPos[0] > xf_[1])
                {
                    analyticSolution_[index] = swr_;
                }
            }

            // non-linear material law

            else
            {
                // find x_f next to global coordinate of the vertex
                int xnext = 0;
                for (int i = intervalNum_; i >= 0; i--)
                {
                    if (globalPos[0] < xf_[i])
                    {
                        xnext = i;
                        break;
                    }
                }

                // account for the area not yet reached by the front
                if (globalPos[0] > xf_[xhelp2])
                {
                    analyticSolution_[index] = swr_;
                    continue;
                }

                if (globalPos[0] <= xf_[xhelp2])
                {
                    analyticSolution_[index] = satVec_[xnext];
                    continue;
                }
            }
        }

        // call function to calculate the saturation error_
        calcSatError(problem_.variables().saturation());

        return;
    }

    void calculateAnalyticSolution()
    {
        time_ += problem_.timeManager().timeStepSize();

        updateExSol();
    }

    //Write saturation and pressure into file
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        BlockVector *analyticSolution = writer.template createField<Scalar, 1> (size_);
        BlockVector *error = writer.template createField<Scalar, 1> (size_);

        *analyticSolution = analyticSolution_;
        *error = error_;

        writer.addCellData(analyticSolution, "saturation (exact solution)");
        writer.addCellData(error, "error_");

        return;
    }

    //! Construct an IMPES object.
    BuckleyLeverettAnalytic(Problem& problem, Scalar totalVelocity = 3e-7) :
        problem_(problem), analyticSolution_(0), error_(0), elementVolume_(0), size_(problem.gridView().size(0)), vTot_(totalVelocity), dummyElement_(
                *(problem_.gridView().template begin<0> ())), dummyGlobal_(GlobalPosition(1))
    {
        initializeAnalytic();
        prepareAnalytic();
    }

protected:


private:
    Problem& problem_;

    BlockVector analyticSolution_;
    BlockVector error_;
    BlockVector elementVolume_;

    Scalar time_;
    int size_;
    Scalar swr_;
    Scalar snr_;
    Scalar porosity_;
    Scalar vTot_;
    enum
    {
        intervalNum_ = 1000, pointNum_ = intervalNum_ + 1
    };
    Dune::FieldVector<Scalar, pointNum_> satVec_;
    Dune::FieldVector<Scalar, pointNum_> fractionalW_;
    Dune::FieldVector<Scalar, pointNum_> dfwdsw_;
    Dune::FieldVector<Scalar, pointNum_> xf_;
    int dfwdswmax_;
    const Element& dummyElement_;
    const GlobalPosition& dummyGlobal_;

};
}
#endif
