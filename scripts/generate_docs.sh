#!/usr/bin/env sh

set -eux

script_dir="$(realpath "$(dirname "$0")")"
cd "${script_dir}"

python3 -m pdoc --docformat numpy --math -o "${script_dir}/../docs/" coffe

cd -

printf "Docs generated at: %s\n" "$(realpath "${script_dir}/../docs/index.html")"

set +eux
