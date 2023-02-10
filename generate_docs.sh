#!/usr/bin/env sh

set -eux

python3 -m pdoc --docformat numpy --math -o docs/ ./coffe

printf "Docs generated at: %s\n" "${PWD}/docs/index.html"

set +eux
