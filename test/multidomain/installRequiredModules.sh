#!/bin/sh

# First create and enter a folder for the dune and dumux modules

# download Dune core modules
git clone http://git.dune-project.org/repositories/dune-common
cd dune-common
git checkout releases/2.3
cd ..

git clone http://git.dune-project.org/repositories/dune-geometry
cd dune-geometry
git checkout releases/2.3
cd ..

git clone http://git.dune-project.org/repositories/dune-grid
cd dune-grid
git checkout releases/2.3
cd ..

git clone http://git.dune-project.org/repositories/dune-istl
cd dune-istl
git checkout releases/2.3
cd ..

git clone http://git.dune-project.org/repositories/dune-localfunctions
cd dune-localfunctions
git checkout releases/2.3
cd ..

# download dune-PDELab
git clone http://git.dune-project.org/repositories/dune-pdelab
cd dune-pdelab
git checkout releases/1.1
cd ..

# download dune-multidomaingrid
git clone git://github.com/smuething/dune-multidomaingrid.git
cd dune-multidomaingrid
git checkout releases/2.3
cd ..

# download dune-multidomain
git clone git://github.com/smuething/dune-multidomain.git
cd dune-multidomain
git checkout deac3cecfc6697c1f5316d55c0fadd74f51d92bc
cd ..

# download DuMuX
svn co svn://svn.iws.uni-stuttgart.de/DUMUX/dumux/trunk dumux

# apply patches to PDELab
cd dune-pdelab
patch -p0 < ../dumux/patches/pdelab-1.1.0.patch
cd ..
