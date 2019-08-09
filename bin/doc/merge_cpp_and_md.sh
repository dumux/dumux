#!/bin/sh

# The arguments should be names of Markdown and/or C++ files.
# The script loops through all arguments. If a Markdown file is encountered,
# its content is forwarded to stdout. For another file, a Markdown header
# stating the filename is printed and the script cpp_to_md.sh is applied.

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
