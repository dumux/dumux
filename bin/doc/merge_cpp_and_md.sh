#!/bin/sh

# check if help is needed
if test "$1" = "--help" || test "$1" = "-help" \
   || test "$1" = "help" || test "$1" = ""; then
  echo ""
  echo "USAGE: $0 FILENAME_1 [FILENAME_2 ...]"
  echo ""
  echo "The arguments should be names of Markdown and/or C++ files."
  echo "The script loops through all arguments. If a Markdown file is encountered,"
  echo "its content is forwarded to stdout. For another file, a Markdown header"
  echo "stating the filename is printed and the script cpp_to_md.sh is applied."
  echo "To be used for creating README.mds in the example folder. For example, by"
  echo "bash $0 intro.md spatialparams.hh problem.hh main.cc results.md >README.md"
  exit 0
fi

for filename in "$@"
do
  if [[ ${filename: -3} == ".md" ]]; then
    cat $filename
  else
    echo ""
    echo ""
    echo "## The file \`$filename\`"
    echo ""
    bash cpp_to_md.sh $filename
    echo ""
  fi
done
