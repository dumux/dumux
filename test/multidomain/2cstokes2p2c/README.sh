#!/bin/sh

# First create and enter a folder for the dune and dumux modules
# If desired, you can enter this folder and execute this script for a DUNE and DUMUX installation.
# For that purpose, you can comment/uncomment required lines below, 
# currently only the DUNE modules are checked out.
# Following dune components should be available in the described versions:

# dune-multidomain (master branch):
git clone git://gitorious.org/dune-multidomain/dune-multidomain.git
cd dune-multidomain
git checkout deac3cecfc6697c1f5316d55c0fadd74f51d92bc
cd ..

# dune-multidomaingrid (master branch):
git clone git://gitorious.org/dune-multidomaingrid/dune-multidomaingrid.git
cd dune-multidomaingrid
git checkout 30ff14d6b49c8adabf9e5cec67f20fcb3270a77e
cd ..

# dune-pdelab (master branch):
git clone http://git.dune-project.org/repositories/dune-pdelab
cd dune-pdelab
git checkout 1de21a0859e10914ba5c3b70f67c5517573ac70e
cd ..

# tested DUNE modules:
svn checkout -r 7436 https://svn.dune-project.org/svn/dune-common/trunk dune-common
svn checkout -r 489 https://svn.dune-project.org/svn/dune-geometry/trunk dune-geometry
svn checkout -r 8930 https://svn.dune-project.org/svn/dune-grid/trunk dune-grid
svn checkout -r 1910 https://svn.dune-project.org/svn/dune-istl/trunk dune-istl
svn checkout -r 1206 https://svn.dune-project.org/svn/dune-localfunctions/trunk dune-localfunctions

# DUMUX
#ssh-add
# for an external installation, use
# svn co --username=anonymous --password=’ ’ svn://svn.iws.uni-stuttgart.de/DUMUX/dumux/trunk dumux

#svn checkout svn+ssh://luftig/home/svn/DUMUX/dumux/trunk dumux-stable
#svn checkout svn+ssh://luftig/home/svn/DUMUX/dune-mux/trunk dumux-devel
#svn checkout svn+ssh://luftig/home/svn/DUMUX/external/trunk external

# patch m4 folder in dumux-stable for compatibility with the DUNE modules above
#patch -p0 < patches/dumux-m4.patch

# install the external modules: UG and SuperLU or PARDISO are required
#cd external
#./installExternal.sh ug
#./installExternal.sh blas
#./installExternal.sh superlu
#cd ..

# dune-pdelab/dune/pdelab/backend/istlvectorbackend.hh has to be replaced by the one in the subfolder additionalfiles
cp dumux-stable/test/multidomain/2cnistokes2p2cni/additionalfiles/istlvectorbackend.hh dune-pdelab/dune/pdelab/backend/

# in dumux/freeflow/new_stokes/stokeslocalresidual:
# if the pressure is set on the entire right boundary, the interpolation at the lower right corner has to be switched off;
# therefore, comment in interpolateCornerPoints the local vertex 1 (lower right corner of Stokes domain);

# run dunecontrol
# ./dune-common/bin/dunecontrol PATH_TO_OPTS all

# the problem hardly converges with BiCGSTAB; use a direct solver like SuperLU or PARDISO
