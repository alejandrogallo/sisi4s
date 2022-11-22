#!/usr/bin/env bash

find src/ -name '*.cxx' -or -name '*.hpp' -exec clang-format -i {} \;
