# One click install script for dumux
echo " "
echo " "
echo "*********************************************************************************************"
echo "(0/3) Checking all prerequistes. (git gcc g++ cmake pkg-config paraview)"
echo "*********************************************************************************************"

# check some prerequistes
for PRGRM in git gcc g++ cmake pkg-config paraview; do
    if ! [ -x "$(command -v $PRGRM)" ]; then
        echo "Error: $PRGRM is not installed." >&2
        exit 1
    fi
done

currentver="$(gcc -dumpversion)"
requiredver="7"
if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" != "$requiredver" ]; then
    echo "gcc greater than or equal to $requiredver is required for dumux releases >=3.2!" >&2
    exit 1
fi

if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "(0/3) An error occured while checking for prerequistes."
    echo "*********************************************************************************************"
    exit $?
else
    echo "*********************************************************************************************"
    echo "(1/3) All prerequistes found."
    echo "*********************************************************************************************"
fi


# make a new folder containing everything
mkdir $(pwd)/DUMUX
cd DUMUX

echo "*********************************************************************************************"
echo "(1/3) Cloning repositories. This may take a while. Make sure to be connected to the internet."
echo "*********************************************************************************************"
DUNE_VERSION=2.7
DUMUX_VERSION=3.2
# the core modules
for MOD in common geometry grid localfunctions istl; do
    git clone -b releases/$DUNE_VERSION  https://gitlab.dune-project.org/core/dune-$MOD.git
done

# dumux
git clone -b releases/$DUMUX_VERSION https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git

if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "(1/3) Failed to clone the repositories. Look for repository specific errors."
    echo "*********************************************************************************************"
    exit $?
else
    echo "*********************************************************************************************"
    echo "(2/3) All repositories have been cloned into a containing folder."
    echo "*********************************************************************************************"
fi

echo " "

echo "**************************************************************************************************"
echo "(2/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes."
echo "**************************************************************************************************"

# run dunecontrol
./dune-common/bin/dunecontrol --opts=cmake.opts all

if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "(2/3) Failed to build the dune libaries."
    echo "*********************************************************************************************"
    exit $?
else
    echo "*****************************************************************************************************"
    echo "(3/3) Succesfully configured and built dune and dumux."
    echo "*****************************************************************************************************"
fi
