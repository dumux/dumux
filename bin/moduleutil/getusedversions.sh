# !/bin/bash

# skript to extract Git version information
# including automating patch extraction and handling of uncommitted changes
#
# (c) 2016 Thomas Fetzer
# (c) 2016 Christoph Grüninger
#

if [ "$1" = "-h" ]; then
  echo "USAGE: ./getDumuxDuneVersions.sh"
  echo; exit
fi

OUTFILE=versionNumbers.txt

# create patches from get local commits
# $1 directory /which is local copy of Git repository to check, assumed to be identical with module name
# $2 script file to append command for checking out Git repository
function getPatchesGit
{
  # if not yet existing, create common patch directory
  if [ ! -d ../patches ]; then
    mkdir ../patches
  fi
  # create directory for patches of module
  mkdir ../patches-$1

  # create patch for local commits
  git format-patch --output-directory ../patches-$1 @{upstream}
  # create extra patch for uncommitted changes if necessary
  if [[ `git status --porcelain | grep -v "build.*"` ]]; then
    git diff > ../patches-$1/9999-uncommitted-changes.patch
  fi

  # move and rename all patches to common patch directory
  MODULE_PWD=$PWD
  cd ../patches-$1 ;
  for file in *.patch
  do
    mv $file ../patches/$1_$file
    printf "Created patch %s\n" $1_$file
    echo "patch -p1 < ../patches/$1_$file" >> $2
  done
  cd $MODULE_PWD
  # remove helper dirctory
  rm -r ../patches-$1
}

# get version from Git, output and write to file
# $1 directory which is local copy of Git repository to check, assumed to be identical with module name
# $2 script file to append command for checking out Git repository
function getVersionGit
{
  # enter directory
  cd $1
  # get URL from remote a.k.a. the repository
  GIT_REMOTE="$(git ls-remote --get-url)"
  # get revision, date, author and branch of last remote/upstream commit
  GIT_UPSTREAM_REVISION="$(git log -n 1 --format=%H @{upstream})"
  GIT_UPSTREAM_DATE="$(git log -n 1 --format=%ai @{upstream})"
  GIT_UPSTREAM_AUTHOR="$(git log -n 1 --format=%an @{upstream})"
  GIT_UPSTREAM_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
  printf "%-25s %-15s %-15s\n" $1 $GIT_UPSTREAM_BRANCH $GIT_UPSTREAM_REVISION
  echo "# $1" >> $2
  echo "# $GIT_UPSTREAM_BRANCH # $GIT_UPSTREAM_REVISION # $GIT_UPSTREAM_DATE # $GIT_UPSTREAM_AUTHOR" >> $2
  echo "git clone $GIT_REMOTE" >> $2
  echo "cd $1" >> $2
  echo "git checkout $GIT_UPSTREAM_BRANCH" >> $2
  echo "git reset --hard $GIT_UPSTREAM_REVISION" >> $2
  # check for untracked files and error out if found, but ignore build*
  if [[ `git ls-files --others --exclude-standard | grep -v "build*" | grep -v "\.[A-Za-z]*"` ]]; then
    echo "Error: Found untracked files in module $1. Please commit, stash or remove them."
    echo "You might want to delete the folder patches/ and the module just created."
    exit
  fi
  # create patches if number of commits is not zero
  if [[ "$(git rev-list HEAD ^@{upstream} --count)" > 0 ]] || [[ `git status --porcelain | grep -v "build.*" | grep -v "\.[A-Za-z]*"` ]]; then
    getPatchesGit $1 $2
  fi
  echo "cd .." >> $2
  echo "" >> $2
  cd ..
}

# run script from command line
# suppressed for use of external script when variable set accordingly
if [ "$CALL_FROM_EXTERNAL_SCRIPT" != "yes" ]; then
  echo "# DUNE/DUMUX VERSIONS" > $OUTFILE

  echo "Creating file containing the version numbers:"
  for dir in dune-* dumux*
  do
    DIR=${dir%*/}
  #   echo $DIR
    # skip folders that are no Git repository
    if [ -d $DIR/.git ]; then
      getVersionGit $DIR ../$OUTFILE
    fi
  done

  echo "done"

  exit
fi
