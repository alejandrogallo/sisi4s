#!/usr/bin/env bash

set -eu

type -f clang-format || {
  echo "install clang-format"
  exit 1
}

fmt="clang-format -i {}"
vendor="-not -path *vendor*"
set -x
find src/ $vendor -name '*.cxx' -exec $fmt \; \
      -or $vendor -name '*.hpp' -exec $fmt \;
