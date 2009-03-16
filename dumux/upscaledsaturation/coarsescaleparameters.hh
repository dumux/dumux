// $Id:$
#ifndef DUNE_COARSESCALEPARAMETERS_HH
#define DUNE_COARSESCALEPARAMETERS_HH

#include<dumux/upscaledsaturation/variablematrix.hh>

namespace Dune
{
template<class Scalar,class DispersionEntryType,class DispersionSatEntryType,class fluxCorrectionEntryType>class CoarseScaleParameters
{
public:
  CoarseScaleParameters()
  {}

      virtual VariableMatrix<Scalar, DispersionEntryType>& getDispersion()
      {
          return dispersion_;
      }

      virtual VariableMatrix<Scalar, DispersionSatEntryType>& getDispersionSat()
      {
          return dispersionSat_;
      }

      virtual VariableMatrix<Scalar, fluxCorrectionEntryType>& getFluxCorr()
      {
          return fluxCorrection_;
      }

      virtual VariableMatrix<Scalar, fluxCorrectionEntryType>& getFluxCorrSat()
      {
          return fluxCorrectionSat_;
      }

private:
    VariableMatrix<Scalar, fluxCorrectionEntryType> fluxCorrection_;
    VariableMatrix<Scalar, fluxCorrectionEntryType> fluxCorrectionSat_;
    VariableMatrix<Scalar, DispersionEntryType> dispersion_;
    VariableMatrix<Scalar, DispersionSatEntryType> dispersionSat_;
};
}
#endif
