### OPM
Within the OPM release 2021.10, some files in opm-common miss include statements. While on some systems such as Ubuntu 20.04, the corresponding headers are included implicitly, explicit includes are required on other systems such as Ubuntu 21.10. To add the includes, apply the patch

 - `opm_common_2021_10.patch` in your directory containing opm-common
   ```bash
   patch -p1 <../dumux/patches/opm_common_2021_10.patch
   ```

### dune-istl
If the AMGBackend should be used without SuperLU as coarse grid solver, it can
be benefitial to decrease the corresponding tolerance. To do so, apply the patch

 - `istl-2.6.patch` in your directory containing dune-istl 2.6
   ```bash
   patch -p1 <../dumux/patches/istl-2.6.patch
   ```
