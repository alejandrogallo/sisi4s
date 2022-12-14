#!/usr/bin/env python3

from testis import read_yaml, compare_energies

out = read_yaml("sisi4s.out.yaml")

compare_energies("correct.out.yaml", "sisi4s.out.yaml", accuracy=1e-7)
