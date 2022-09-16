#!/usr/bin/env bash

set -eu
cd `git rev-parse --show-toplevel`

type -a autoreconf > /dev/null ||
{
        cat <<EOF && exit

You don't seem to have autotools installed, please install it.

  - https://www.gnu.org/software/autoconf/
  - https://www.gnu.org/software/automake/

EOF
}


cat <<EOF

  Creating configure script

EOF


autoreconf -vif .
test -f configure || {
  cat <<EOF

  An error happened and a configure script could not be built!

EOF
  exit 1
}


cat <<EOF

  Now you can build by doing

  mkdir build
  cd build
  ../configure
  make extern
  make all

EOF
