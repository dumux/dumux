# Benchmark Richards equation evaportion

The benchmark case is a benchmarking project for hydrological models @cite Vanderborght2005.
It was also part of the collaborative benchmark project on
root-soil interaction [2] in which DuMux participated as a simulator.

## Benchmark and model description

The Richards equation is used to model the unsaturated flow in the soil.


We solve the Richards equation,


## Reference solution

The reference solution is an approximate analytical solution.

The evaporation rate is given by

Evaporationrate in mm/day (Eq. 46 & 47 in @cite Vanderborght2005)


## Reporting


## How to reproduce

```sh
cd <build-dir>/test/porousmediumflow/richards/benchmarks
make test_richards_benchmark_evaporation_tpfa
python run_and_plot_m22.py
```

## Result

![evaporation_benchmark_results](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/test/porousmediumflow/richards/benchmarks/doc/evaporation_benchmark_results.svg)
