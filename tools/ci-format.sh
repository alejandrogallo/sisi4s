#!/usr/bin/env bash

set -eu

type -f clang-format || {
  echo "install clang-format"
  exit 1
}

fmt="clang-format -i {}"
find src/ -name '*.cxx' -exec $fmt \; \
      -or -name '*.hpp' -exec $fmt \;
