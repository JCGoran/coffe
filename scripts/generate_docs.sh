#!/usr/bin/env sh

set -eux

# dir where the script is
script_dir="$(cd "$(dirname "$0")"; pwd -P)"

python -m pdoc --docformat numpy --math -o "${script_dir}/../docs/" coffe

printf "Docs generated at: %s\n" "${script_dir}/../docs/index.html"

set +eux
