#include <algorithms/XYZReader.hpp>
#include <util/Libint.hpp>
#include <util/Log.hpp>

#define LOGGER(_l) LOG(_l, "XyzReader")

using namespace sisi4s;

using Structure = std::vector<libint2::Atom>;

DEFSPEC(
    XyzReader,
    SPEC_IN({"path", SPEC_VALUE_DEF("Path to a xyz file", std::string, "")},
            {"structure",
             SPEC_VALUE_DEF("Xyz File directly as a string", std::string, "")}),
    SPEC_OUT({"atoms",
              SPEC_VAROUT("A vector of libint atoms defining a Xyz structure",
                          Structure *)}));

IMPLEMENT_EMPTY_DRYRUN(XyzReader) {}

IMPLEMENT_ALGORITHM(XyzReader) {

  const std::string xyzStructureFile = in.get<std::string>("path"),
                    structure_string = in.get<std::string>("structure");

  Structure atoms;
  if (xyzStructureFile.size()) {
    LOGGER(1) << "structure: " << xyzStructureFile << std::endl;
    std::ifstream structureFileStream(xyzStructureFile.c_str());
    atoms = libint2::read_dotxyz(structureFileStream);
    structureFileStream.close();
  } else if (structure_string.size()) {
    std::istringstream s;
    s.str(structure_string);
    atoms = libint2::read_dotxyz(s);
  } else {
    throw "path or structure has to be provided";
  }

  LOGGER(1) << "#atoms: " << atoms.size() << std::endl;

  out.set<Structure *>("atoms", new Structure(atoms));
}
