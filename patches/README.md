### OPM
OPM release 2024.10 is incompatible with Dune > 2.9 since deprecated names have been removed there. To switch to the new names, apply the patch

 - `opm_grid_2024_10.patch` in your directory containing opm-grid
   ```bash
   patch -p1 <../dumux/patches/opm_grid_2024_10.patch
   ```

### dune-istl
If the AMGBackend should be used without SuperLU as coarse grid solver, it can
be benefitial to decrease the corresponding tolerance. To do so, apply the patch

 - `istl-2.6.patch` in your directory containing dune-istl 2.6
   ```bash
   patch -p1 <../dumux/patches/istl-2.6.patch
   ```
