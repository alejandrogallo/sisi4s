#include <string>
#include <algorithms/NuclearRepulsionEnergy.hpp>
#include <util/Libint.hpp>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <util/Emitter.hpp>

using namespace sisi4s;

DEFSPEC(NuclearRepulsionEnergy,
        SPEC_IN({"atoms",
                 SPEC_VARIN(
                     "Vector of libint atoms specifying a molecular structure",
                     std::vector<libint2::Atom> *)}),
        SPEC_OUT());

IMPLEMENT_EMPTY_DRYRUN(NuclearRepulsionEnergy) {}
#define LOGGER(_l) LOG(_l, "NuclearRepulsionEnergy")

static double getEnergy(const std::vector<libint2::Atom> &structure) {
  unsigned int i, j;
  double enuc(0.0), r2(0.0);
  LOGGER(1) << "Calculating nuclear repulsion energy" << std::endl;
  for (i = 0; i < structure.size(); ++i) {
    for (j = i + 1; j < structure.size(); ++j) {
      r2 = 0.0;
      r2 += pow(structure[i].x - structure[j].x, 2);
      r2 += pow(structure[i].y - structure[j].y, 2);
      r2 += pow(structure[i].z - structure[j].z, 2);
      enuc +=
          structure[i].atomic_number * structure[j].atomic_number / sqrt(r2);
    }
  }
  return enuc;
}

IMPLEMENT_ALGORITHM(NuclearRepulsionEnergy) {

  auto const &atoms = *in.get<std::vector<libint2::Atom> *>("atoms");

  LOGGER(1) << "#atoms: " << atoms.size() << std::endl;

  double enuc = getEnergy(atoms);

  LOGGER(1) << std::setprecision(15) << std::setw(10) << "energy = " << enuc
            << std::endl;

  EMIT() << YAML::Key << "energy" << YAML::Value << enuc;
}
