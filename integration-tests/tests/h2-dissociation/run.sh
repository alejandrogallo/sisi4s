lengths=(
  0.00001
  0.0001
  0.001
  0.006
  0.011
  0.016
  0.021
  0.026
  0.031
  0.036
  0.041
  0.046
  0.051
  0.056
  0.061
  0.066
  0.071
  0.076
  0.081
  0.086
  0.091
  0.096
  0.1
  0.2
  0.3
  0.4
  0.5
  0.7
  0.8
  1.0
  1.2
  1.4
  1.6
  1.8
  2.0
  2.2
  2.4
  2.6
  2.8
  3.0
  3.2
  3.4
  3.6
  3.8
  4.0
  4.2
  4.4
  4.6
  4.8
  5.0
)

for length in ${lengths[@]}; do
  m4 -D MOL_LENGTH=${length} ./in.yaml.template >in."${length}".yaml
  mpirun \
    ./sisi4s \
    -i in."${length}".yaml \
    -o out."${length}".yaml |
      tee out."${length}".stdout
done
