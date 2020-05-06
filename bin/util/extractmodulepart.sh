#! /bin/bash

# check if help is needed
if test "$1" = "--help" || test "$1" = "-help" \
   || test "$1" = "help" || test "$1" = ""; then
  echo ""
  echo "USAGE: $0 MODULE_DIR FOLDER_1 [FOLDER_2 ...]"
  echo ""
  echo "MODULE_DIR is the folder containing the DUNE module from which you"
  echo "want to extract. The script has to be called one level above it."
  echo ""
  echo "The FOLDERs need to indicate subfolders of MODULE_DIR. At least one"
  echo "of them has to contain a source file *.cc of an executable for which"
  echo "you would like to timber a table in dumux-pub."
  exit 0
fi

MODULE_DIR=$1
# if MODULE_DIR contains a slash as last character, delete it
MODULE_DIR=$(echo $MODULE_DIR | sed s'/\/$//')
# check if we are above MODULE_DIR
if (! [ -d $MODULE_DIR ]); then
  echo "ERROR: you need to run the script one level above the folder $MODULE_DIR."
  echo "Run \"$0 --help\" for details."
  exit 1
fi

# determine all source files in the paths passed as arguments
cd $MODULE_DIR
MODULE_FULL_PATH=$(pwd)
ALL_SOURCES=""
ALL_DIRECTORIES=""
for DIR_PATH in ${*:2}; do
  STRIPPED_PATH=$(sed 's%.*'"$MODULE_DIR"'/%%g' <(echo "$DIR_PATH"))
  ALL_DIRECTORIES="$ALL_DIRECTORIES $STRIPPED_PATH"
  ALL_SOURCES="$ALL_SOURCES $(find $STRIPPED_PATH -name '*.cc' 2>/dev/null)"
done
cd ..

# check if sources have been obtained
CONTRACTED="$(echo "$ALL_SOURCES" | tr -d " \t\n\r")"
if test "$CONTRACTED" = ""; then
  echo "ERROR: no source files *.cc found in the directories ${*:2}."
  echo "Be sure to provide a list of paths as arguments to this script."
  echo "Run \"$0 --help\" for details."
  exit 1;
fi

# try to find the duneproject script
if hash duneproject 2>/dev/null; then
  DUNE_PROJECT="duneproject"
else
  DUNE_PROJECT=$(find . -name 'duneproject')
fi
if test "$DUNE_PROJECT" = ""; then
  echo "ERROR: Could not find duneproject."
  echo "Be sure to either have duneproject in your search path"
  echo "or to run this script from a directory that contains duneproject."
  exit 1;
fi

# give explanations
echo ""
echo "This script will"
echo "- extract the following sub-folders of $MODULE_DIR:"
echo ""
for DIR_PATH in $ALL_DIRECTORIES; do
  echo "  $DIR_PATH,"
done
echo ""
echo "  and all headers in $MODULE_DIR that are required to build the"
echo "  executables from the sources"
echo ""
for SOURCE in $ALL_SOURCES; do
  echo "  $SOURCE,"
done
echo ""
echo "- copy the extracted files into a freshly created DUNE module, retaining the"
echo "  directory structure,"
echo ""
echo "- update/create all required files like CMakeLists.txt,"
echo ""
echo "- store the versions of all used Dune module"
echo ""
echo "- and extract their modifications as patches."
echo ""
echo "Thus, you receive a fully-working DUNE module containing the subset of"
echo "$MODULE_DIR that is required to run your application."
echo ""
echo "duneproject will be run now. The new module should NOT depend on the"
echo "module in $MODULE_DIR."
echo ""
read -p "Read the above and press [Enter] to proceed..."

# run duneproject
OLD_LS="$(ls)"
$DUNE_PROJECT
NEW_LS="$(ls)"

# determine the new module/directory name
DIFF_OUTPUT=$(diff <(echo "$OLD_LS" ) <(echo "$NEW_LS"))
FOUND="0"
MODULE_NAME=""
for WORD in $DIFF_OUTPUT; do
  if test "$FOUND" = "1"; then
    MODULE_NAME=$WORD
  fi
  if test "$WORD" = ">"; then
    FOUND="1"
  fi
done
if test "$MODULE_NAME" = ""; then
  echo "ERROR: could not find new module. Aborting."
  exit 1
else
  echo ""
  echo "$(basename $0): Found new module $MODULE_NAME"
fi

echo -n "Determining required headers..."
cd $MODULE_DIR

# initialize lists to hold required headers
LAST_REQUIRED_HEADERS=""
REQUIRED_HEADERS="$ALL_SOURCES"

while test "$LAST_REQUIRED_HEADERS" != "$REQUIRED_HEADERS"; do
  echo -n "."
  LAST_REQUIRED_HEADERS="$REQUIRED_HEADERS"

  # remove the file that stores all required headers line by line
  rm -f tmp_header_file

  for HEADER in $REQUIRED_HEADERS; do
    # if tmp_header_file does not exist, create it
    if test "$(ls tmp_header_file 2>/dev/null)" = ""; then
      touch tmp_header_file
    fi

    # append header name to tmp_header_file
    echo "$HEADER" >> tmp_header_file

    # obtain the base name and the directory name of the header that is considered
    HEADER_BASE_NAME=$(basename $HEADER)
    HEADER_DIR_NAME=$(dirname $HEADER)

    # create a list of all files that are included from the header
    sed -i 's/#include/#include /' $HEADER
    INCLUDE_LIST=$(tr -s " " < $HEADER | tr -d '><"' | awk '/#include/{print $2}')
    sed -i 's/#include /#include/' $HEADER

    # look at each of the included files
    for INCLUDED_HEADER in $INCLUDE_LIST; do
      # don't include config.h
      if test "$INCLUDED_HEADER" = "config.h"; then
        continue
      fi

      INCLUDED_BASE_NAME=$(basename $INCLUDED_HEADER)
      INCLUDED_DIR_NAME=$(dirname $INCLUDED_HEADER)

      # if the header file exists, add it
      if test "$(ls $INCLUDED_HEADER 2>/dev/null)" = "$INCLUDED_HEADER"; then
        echo "$INCLUDED_HEADER" >> tmp_header_file
        continue
      fi

      # try the header preceded by its path
      INCLUDED_HEADER_WITH_PATH="${HEADER_DIR_NAME}/${INCLUDED_HEADER}"

      # if the header file actually exists, add it
      if test "$(ls $INCLUDED_HEADER_WITH_PATH 2>/dev/null)" = "$INCLUDED_HEADER_WITH_PATH"; then
        # remove "../" from the path
        cd $(dirname $INCLUDED_HEADER_WITH_PATH)
        HEADER_FULL_PATH=$(pwd)
        HEADER_RELATIVE_PATH=${HEADER_FULL_PATH#$MODULE_FULL_PATH}
        HEADER_RELATIVE_PATH=$(echo $HEADER_RELATIVE_PATH | sed 's/^.//')
        INCLUDED_HEADER_WITH_PATH="${HEADER_RELATIVE_PATH}/${INCLUDED_BASE_NAME}"
        cd $MODULE_FULL_PATH
        echo "$INCLUDED_HEADER_WITH_PATH" >> tmp_header_file
      fi
    done
  done

  # sort the required headers, eliminate copies
  REQUIRED_HEADERS=$(sort -u tmp_header_file)

done

# remove the file that stores all required headers line by line
rm -f tmp_header_file

# provide some output and copy everything to the new module
echo ""
echo -n "Number of required headers: "
echo "$REQUIRED_HEADERS" | wc -w
for HEADER in $REQUIRED_HEADERS; do
  echo $HEADER
  rsync -R $HEADER ../$MODULE_NAME
done
for DIR_PATH in $ALL_DIRECTORIES; do
  rsync -rR $DIR_PATH ../$MODULE_NAME
done

# delete all architecture-dependent files
cd ../$MODULE_NAME
find . -name Makefile.in -delete
find . -name Makefile -delete
find . -name '*.o' -delete
find . -wholename '*.deps/*' -delete
find . -path '*.deps' -delete
echo "Removed architecture-dependent files."

# remove directories that are not required
sed -i '/(dune)/d' CMakeLists.txt
sed -i '/(src)/d' CMakeLists.txt
rm -rf dune/
rm -rf src/

# create a list of the subfolders
EMPTY_DIR_NAME="."
rm -f tmp_subdir_file
for HEADER in $REQUIRED_HEADERS; do

  # move through every header, cutting off the last folder
  while test "$HEADER" != $EMPTY_DIR_NAME; do
    if test "$(dirname $HEADER)" != $EMPTY_DIR_NAME; then
      dirname $HEADER >> tmp_subdir_file
    fi

    HEADER=$(dirname $HEADER)
  done

  # finally add the topmost folder
  if test "$(basename $HEADER)" != $EMPTY_DIR_NAME; then
    basename $HEADER >> tmp_subdir_file
  fi
done
SUBDIR_LIST=$(sort -u tmp_subdir_file)
rm -f tmp_subdir_file

# treat every subfolder
START_DIR=$PWD
INITIAL_SUBDIR_LIST=""
for SUBDIR in $SUBDIR_LIST; do
  if test "$SUBDIR" = "$(basename $SUBDIR)"; then
    INITIAL_SUBDIR_LIST="$INITIAL_SUBDIR_LIST $SUBDIR"
  fi

  cd $SUBDIR

  # create lists of files and folders
  FILE_LIST=""
  DIR_LIST=""
  for FILE_OR_FOLDER in $( ls ); do
    if ([ -f $FILE_OR_FOLDER ]); then
      FILE_LIST="$FILE_LIST $FILE_OR_FOLDER"
    fi

    if ([ -d $FILE_OR_FOLDER ]); then
      DIR_LIST="$DIR_LIST $FILE_OR_FOLDER"
    fi
  done

  # create CMakeLists.txt if not present
  if ([ ! -f CMakeLists.txt ]); then
    # add subfolder information
    for DIR in $DIR_LIST; do
      echo "add_subdirectory($DIR)" >> CMakeLists.txt
    done

    # determine whether to add file information
    ADD_FILE_INFO="0"
    for FILE in $FILE_LIST; do
      if ([ -f $FILE ]); then
        ADD_FILE_INFO="1"
      fi
    done

    # add file information
    if test "$ADD_FILE_INFO" = "1"; then
      echo "" >> CMakeLists.txt
      echo "install(FILES"  >> CMakeLists.txt

      for FILE in $FILE_LIST; do
        echo "        $FILE" >> CMakeLists.txt
      done

      echo "        DESTINATION \${CMAKE_INSTALL_INCLUDEDIR}/$SUBDIR)" >> CMakeLists.txt
    fi
  fi

  cd $START_DIR

done

# update top CMakeLists.txt
for INITIAL_SUBDIR in $INITIAL_SUBDIR_LIST; do
  sed -i '/add_subdirectory(doc)/a add_subdirectory('$INITIAL_SUBDIR')' CMakeLists.txt
done
echo "Updated and created CMakeLists.txt's."
cd ..

# add includes to the automatically generated file *Macros.cmake
MACROS_FILE="$(find $MODULE_NAME -name '*Macros.cmake')"
if test "$(find $MODULE_NAME -type f -name 'CMakeLists.txt' -exec grep add_csv_file_links {} \; 2>/dev/null)" != ""; then
   cp dumux-devel/cmake/modules/AddCSVFileLinks.cmake $MODULE_NAME/cmake/modules/.
   echo "include(AddCSVFileLinks)" >>$MACROS_FILE
fi
if test "$(find $MODULE_NAME -type f -name 'CMakeLists.txt' -exec grep add_executable_all {} \; 2>/dev/null)" != ""; then
   cp dumux-devel/cmake/modules/AddExecutableAll.cmake $MODULE_NAME/cmake/modules/.
   echo "include(AddExecutableAll)" >>$MACROS_FILE
fi
if test "$(find $MODULE_NAME -type f -name 'CMakeLists.txt' -exec grep add_file_link {} \; 2>/dev/null)" != ""; then
   cp dumux-devel/cmake/modules/AddFileLink.cmake $MODULE_NAME/cmake/modules/.
   echo "include(AddFileLink)" >>$MACROS_FILE
fi
if test "$(find $MODULE_NAME -type f -name 'CMakeLists.txt' -exec grep add_folder_link {} \; 2>/dev/null)" != ""; then
   cp dumux-devel/cmake/modules/AddFolderLink.cmake $MODULE_NAME/cmake/modules/.
   echo "include(AddFolderLink)" >>$MACROS_FILE
fi
if test "$(find $MODULE_NAME -type f -name 'CMakeLists.txt' -exec grep add_gnuplot_file_links {} \; 2>/dev/null)" != ""; then
   cp dumux-devel/cmake/modules/AddGnuplotFileLinks.cmake $MODULE_NAME/cmake/modules/.
   echo "include(AddGnuplotFileLinks)" >>$MACROS_FILE
fi

# create $README_FILE
README_FILE="README.md"
mv $MODULE_NAME/README $MODULE_NAME/$README_FILE
echo "This file has been created automatically. Please adapt it to your needs." >$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE
echo "===============================" >>$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE
echo "## Content" >>$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE
echo "The content of this DUNE module was extracted from the module \`$MODULE_DIR\`." >>$MODULE_NAME/$README_FILE
echo "In particular, the following subfolders of \`$MODULE_DIR\` have been extracted:" >>$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE
for DIR_PATH in $ALL_DIRECTORIES; do
  echo "*  \`$DIR_PATH\`," >>$MODULE_NAME/$README_FILE
done
echo "" >>$MODULE_NAME/$README_FILE
echo "Additionally, all headers in \`$MODULE_DIR\` that are required to build the" >>$MODULE_NAME/$README_FILE
echo "executables from the sources" >>$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE
for SOURCE in $ALL_SOURCES; do
  echo "*  \`$SOURCE\`," >>$MODULE_NAME/$README_FILE
done
echo "" >>$MODULE_NAME/$README_FILE
echo "have been extracted. You can configure the module just like any other DUNE" >>$MODULE_NAME/$README_FILE
echo "module by using \`dunecontrol\`. For building and running the executables," >>$MODULE_NAME/$README_FILE
echo "please go to the build folders corresponding to the sources listed above." >>$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE

# get versions from modules we are depeding on
# create patches for un-pushed commits and uncomitted changes

# try to find the dunecontrol script
if hash dunecontrol 2>/dev/null; then
  DUNE_CONTROL="dunecontrol"
else
  DUNE_CONTROL=$(find . -type f -executable -name 'dunecontrol')
fi
if test "$DUNE_CONTROL" = ""; then
  echo "ERROR: Could not find dunecontrol."
  echo "Be sure to either have dunecontrol in your search path"
  echo "or to run this script from a directory that contains dunecontrol."
  exit
fi
# get names of all modules we are depeding on
echo "Extracting dependencies to other Dune modules"
DEPENDING_MODULE_NAMES=$($DUNE_CONTROL --module=$MODULE_NAME 2>/dev/null | head -n 1)
DEPENDING_MODULE_NAMES=$(sed 's/--- going to build //' <<< "$DEPENDING_MODULE_NAMES")
DEPENDING_MODULE_NAMES=$(sed 's/  ---//' <<< "$DEPENDING_MODULE_NAMES")
# load script
CALL_FROM_EXTERNAL_SCRIPT="yes"
SCRIPT_DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$SCRIPT_DIR" ]]; then
  SCRIPT_DIR="$PWD";
fi
# create install script
touch $MODULE_NAME/install$MODULE_NAME.sh
. "$SCRIPT_DIR/getusedversions.sh"
# execute function from script
echo "Extracting Git status of dependencies and creating patches"
echo "## Dependencies on other DUNE modules" >>$MODULE_NAME/$README_FILE
echo "" >>$MODULE_NAME/$README_FILE
echo "| module | branch | commit hash |" >> $MODULE_NAME/$README_FILE
echo "|:-------|:-------|:------------|" >> $MODULE_NAME/$README_FILE
for MOD in $DEPENDING_MODULE_NAMES
do
  if [ $MOD != $MODULE_NAME ]; then
    MODULE_DIR=${MOD%*/}
    VERSION_OUTPUT=$(getVersionGit $MODULE_DIR $(pwd)/$MODULE_NAME/install$MODULE_NAME.sh)
    echo "$VERSION_OUTPUT"
    grep "Error:" <<< $VERSION_OUTPUT > /dev/null
    EXIT_CODE=$?
    if [[ $EXIT_CODE == 0 ]]; then
      exit
    fi
    VERSION_OUTPUT=$(echo "$VERSION_OUTPUT" | tr -s " ")
    DEPENDENCY_NAME=$(cut -d " " -f1 <<< $VERSION_OUTPUT)
    BRANCH_NAME=$(cut -d " " -f2 <<< $VERSION_OUTPUT)
    COMMIT_HASH=$(cut -d " " -f3 <<< $VERSION_OUTPUT)
    echo "| $DEPENDENCY_NAME | $BRANCH_NAME | $COMMIT_HASH |" >> $MODULE_NAME/$README_FILE
  fi
done

# move patches folder into module if existing
if [[ -d patches ]]; then
  mv patches $MODULE_NAME
fi

# output guidence for users
echo ""
echo "*************************************************************************"
echo "The extracted module is contained in the subfolder \"$MODULE_NAME\"."
echo "You can build it using \"dunecontrol ... --only=$MODULE_NAME all\"."
echo "*************************************************************************"
echo "BEFORE building, you can add the module to dumux-pub by something like:"
echo "(Rename module name if it does not match the AuthorLastNameYearx scheme"
echo "and commit it to the Git repository dumux-pub)"
echo "git clone https://git.iws.uni-stuttgart.de/dumux-pub/AuthorLastNameYearx.git"
echo "sed -i '/Module:/c\Module: AuthorLastNameYearx' $MODULE_NAME/dune.module"
echo "mv $MODULE_NAME/* dumux-pub/AuthorLastNameYearx/."
echo "cd AuthorLastNameYearx"
echo "git commit -a"
echo "git push"

exit 0
