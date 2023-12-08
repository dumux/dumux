#!/bin/bash
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

options="--error-exitcode=1 --inline-suppr -iformat/fmt/ --enable=missingInclude --suppress=*:*/dune/* --suppress=*:*/fmt/* --suppress=*:*/xml/* --suppress=*:*/json/* --suppress=missingIncludeSystem:*/dumux/* --suppress=preprocessorErrorDirective:*/pyconfig.h"
tmp_dir=$(mktemp -d)
if [[ -s "affectedtests.json" ]]; then
    affected_targets="$(jq -r 'to_entries[] | "\(.value | .target)"' affectedtests.json)"
fi
jq -c '.[]' build-cmake/compile_commands.json | while read -r compile_command; do
    target=$(echo "$compile_command" | sed "s/.*CMakeFiles\/\(.*\)\.dir.*/\1/")
    if ([[ -n $affected_targets ]] && [[ ! $affected_targets =~ $target ]]); then
        continue
    fi
    tmp_file=$(mktemp "${tmp_dir}/XXX.json")
    echo "[$compile_command]" > "$tmp_file"
    echo "$tmp_file"
done | parallel --keep-order cppcheck "$options" --project={}
return_value=$?
rm -r "$tmp_dir"
exit $return_value
