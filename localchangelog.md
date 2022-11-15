- add constantcementmodel for the calculation of mechanical properties of two-component matrix (Rock solid + cement) at different porosities.
- `EffectiveStressLaw` takes one more template class `StressDropLaw`
  `template<class StressType, class StressDropLaw, class GridGeometry>`
  set `StressDropLaw = void` to skip the stress drop calculation
- unit tests in folder `geomechanics/mytests`, `test.py` for automated test
- add new class `StressDropLaw` and `StressDropLawParams` (example implemented in geo/1p test''')
- add MohrSpace  `<geomechanics/math>`, calculations for points, lines and circles.