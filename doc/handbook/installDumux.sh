# One click install script for dumux

# make a new folder containing everything
mkdir $(pwd)/DUMUX
cd DUMUX

echo "*************************************************"
echo "(1/2) Cloning repositories. This may take a while.
Make sure to be connected to the internet."
echo "*************************************************"
# the core modules
for MOD in common geometry grid localfunctions istl; do
    if [ ! -d "dune-$MOD" ]; then
        git clone -b releases/2.6 https://gitlab.dune-project.org/core/dune-$MOD.git
    else
        echo "Skip cloning dune-$MOD because the folder already exists."
        cd dune-$MOD
        git checkout releases/2.6
        cd ..
    fi
done

# dumux
if [ ! -d "dumux" ]; then
    git clone -b releases/3.0 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git
else
    echo "Skip cloning dumux because the folder already exists."
    cd dumux
    git checkout releases/3.0
    cd ..
fi

if [ $? -ne 0 ]; then
    echo "*************************************************"
    echo "Failed to clone the repositories."
    echo "*************************************************"
    exit $?
fi

echo "*************************************************"
echo "(2/2) Configure dune modules and dumux. Build the
dune libaries. This may take several minutes."
echo "*************************************************"
# run build
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all
#
if [ $? -ne 0 ]; then
    echo "*************************************************"
    echo "Failed to build the dune libaries."
    echo "*************************************************"
    exit $?
fi

# echo result
echo "*************************************************"
echo "Successfully configured and built dune and dumux."
echo "Please run the test_dumux.sh script to confirm everything works."
echo "*************************************************"
