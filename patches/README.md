If the AMGBackend should be used without SuperLU as coarse grid solver, it can
be benefitial to decrease the corresponding tolerance. To do so, apply the patch

 - `istl-2.6.patch` in your directory containing dune-istl 2.6
   ```bash
   patch -p1 <../dumux/patches/istl-2.6.patch
   ```
