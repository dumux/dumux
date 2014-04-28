#!/bin/sh

# First create and enter a folder for the dune and dumux modules

# download Dune core modules
git clone http://git.dune-project.org/repositories/dune-common
cd dune-common
git checkout 6fb2492cca04e07ad43074a29667b633b4fa0680
cd ..

git clone http://git.dune-project.org/repositories/dune-geometry
cd dune-geometry
git checkout cda1d514d79f13e70de2d55fdf6906864c2fdcdd
cd ..

git clone http://git.dune-project.org/repositories/dune-grid
cd dune-grid
git checkout cfec4c46bd59219337660ed833fe38cd2acbd364
cd ..

git clone http://git.dune-project.org/repositories/dune-istl
cd dune-istl
git checkout b2c9640a4873ca905caf8e29ebc348a817801e9b
cd ..

git clone http://git.dune-project.org/repositories/dune-localfunctions
cd dune-localfunctions
git checkout d7776478e0e56beceebc48ddfca068b1c115dff3
cd ..

# download dune-PDELab
git clone http://git.dune-project.org/repositories/dune-pdelab
cd dune-pdelab
git checkout releases/1.1
cd ..

# download dune-multidomaingrid
git clone https://users.dune-project.org/repositories/projects/dune-multidomaingrid.git
cd dune-multidomaingrid
git checkout 30ff14d6b49c8adabf9e5cec67f20fcb3270a77e
cd ..

# download dune-multidomain
git clone https://users.dune-project.org/repositories/projects/dune-multidomain.git
cd dune-multidomain
git checkout deac3cecfc6697c1f5316d55c0fadd74f51d92bc
cd ..

# download DuMuX
svn co svn://svn.iws.uni-stuttgart.de/DUMUX/dumux/branches/release-2.5 dumux

# apply patches to PDELab
cd dune-pdelab
patch -p1 < ../dumux/patches/pdelab-1.1.0.patch
cd ..
