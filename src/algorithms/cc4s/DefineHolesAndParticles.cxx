#include <Sisi4s.hpp>
#include <util/Yaml.hpp>
#include <util/Tensor.hpp>
#include <util/Log.hpp>
#include <algorithms/cc4s/DefineHolesAndParticles.hpp>

namespace sisi4s {

IMPLEMENT_EMPTY_DRYRUN(DefineHolesAndParticles) {}

DEFSPEC(DefineHolesAndParticles,
        SPEC_IN({"fileName",
                 SPEC_VALUE("The yaml file path", std::string)->require()},
                {"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},
                {"OrbitalCoefficients",
                 SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({"HoleEigenEnergies",
                  SPEC_VAROUT("The hole eigen energies", Tensor<double> *)},
                 {"ParticleEigenEnergies",
                  SPEC_VAROUT("The hole eigen energies", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(DefineHolesAndParticles) {

  const std::string fileName = in.get<std::string>("fileName");

  cc4s::HPHeader h = YAML::LoadFile(fileName).as<cc4s::HPHeader>();

  LOG(0, "DefineHolesAndParticles") << "fermiEnergy: " << h.fermiEnergy << "\n";

  LOG(0, "DefineHolesAndParticles")
      << "read elements: " << h.energies.size() << "\n";

  const auto it = std::partition(h.energies.begin(),
                                 h.energies.end(),
                                 [&h](double e) { return e <= h.fermiEnergy; });

  const TensorIndex No = it - h.energies.begin(), Nv = h.energies.size() - No,
                    ns = NS;

  auto *epsi = new Tensor<double>(1, &ns, &No),
       *epsa = new Tensor<double>(1, &ns, &Nv);

  LOG(0, "DefineHolesAndParticles") << "No: " << No << "\n";
  LOG(0, "DefineHolesAndParticles") << "Nv: " << Nv << "\n";

  {
    // write epsi
    auto &t = epsi;
    const auto count = Sisi4s::world->rank == 0 ? No : 0;
    std::vector<TensorIndex> indices(count);
    std::iota(indices.begin(), indices.end(), 0);
    t->write(count, indices.data(), h.energies.data());
    out.set<Tensor<double> *>("HoleEigenEnergies", t);
  }

  {
    // write epsa
    auto &t = epsa;
    const auto count = Sisi4s::world->rank == 0 ? Nv : 0;
    std::vector<TensorIndex> indices(count);
    std::iota(indices.begin(), indices.end(), 0);
    t->write(count, indices.data(), h.energies.data() + No);
    out.set<Tensor<double> *>("ParticleEigenEnergies", t);
  }
}

} // namespace sisi4s

namespace YAML {

using namespace sisi4s::cc4s;
template <>
struct convert<HPHeader> {
  using Native = HPHeader;
  static bool decode(Node const &n, Native &r) {
    if (!n.IsMap()) return false;
    YAML_ASSERT_KEY(n, "metaData");
    YAML_ASSERT_KEY(n["metaData"], "fermiEnergy");
    YAML_ASSERT_KEY(n["metaData"], "energies");
    r.fermiEnergy = n["metaData"]["fermiEnergy"].as<double>();
    r.energies = n["metaData"]["energies"].as<std::vector<double>>();
    return true;
  }
};

} // namespace YAML
