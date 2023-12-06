#!/bin/bash
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

options="--error-exitcode=1 --inline-suppr -iformat/fmt/ --enable=missingInclude --suppress=*:*/dune/* --suppress=*:*/fmt/* --suppress=*:*/xml/* --suppress=*:*/json/* --suppress=missingIncludeSystem:*/dumux/*"
tmp_dir=$(mktemp -d)
jq -c '.[]' build-cmake/compile_commands.json | while read -r compile_command; do
    tmp_file=$(mktemp "${tmp_dir}/XXX.json")
    echo "[$compile_command]" > "$tmp_file"
    echo "$tmp_file"
done | parallel --keep-order cppcheck "$options" --project={}
return_value=$?
rm -r "$tmp_dir"
exit $return_value
