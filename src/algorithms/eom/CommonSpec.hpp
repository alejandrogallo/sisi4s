#ifndef EOM_COMMONSPEC_HPP
#define EOM_COMMONSPEC_HPP

#define TODO_EOM_F double
#define EOM_COMMON_SPEC_IN                                                     \
  {"amplitudesConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},        \
      {"type",                                                                 \
       SPEC_ONE_OF("Type of EOM Theory", std::string, "ea", "ip", "ee")},      \
      {"restricted", SPEC_VALUE_DEF("Use restricted equations", bool, false)}, \
      {"energyConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},        \
      {"eigenstates", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},                \
      {"intermediates", SPEC_VALUE_DEF("TODO: DOC", bool, true)},              \
      {"maxIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 32)},             \
      {"minIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},              \
      {"refreshOnMaxBasisSize", SPEC_VALUE_DEF("TODO: DOC", bool, false)},     \
      {"maxBasisSize", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},              \
      {"oneBodyRdmRange", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},       \
      {"printEigenvectorsRange",                                               \
       SPEC_VALUE_DEF("TODO: DOC", std::string, "")},                          \
      {"refreshIterations", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},     \
      {"HHFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},             \
      {"HoleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},        \
      {"HPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},             \
      {"ParticleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},    \
      {"PPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},             \
      {"DoublesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)},    \
      {"HHHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"HHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"HHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"HPHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"HPHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"HPPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"HPPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"PHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
      {"PHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)}, \
  {                                                                            \
    "SinglesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<TODO_EOM_F> *)         \
  }

#endif
