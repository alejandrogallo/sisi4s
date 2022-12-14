#ifndef EMITTER_DEFINED
#define EMITTER_DEFINED

#include <util/SharedPointer.hpp>

#include <string>
#include <yaml-cpp/yaml.h>

namespace sisi4s {
/**
 * \brief Class with static members offering control over yaml emitting.
 * Entries are emitted with the macro EMIT.
 */
class Emitter {
public:
  static void setRank(const int rank);
  static int getRank();
  static YAML::Emitter &getEmitter();
  static void flush();
  static void setFileName(const std::string &);

protected:
  static int rank;
  static std::string fileName;
  static PTR(std::ofstream) yamlFile;
  static PTR(YAML::Emitter) yamlEmitter;
};
} // namespace sisi4s

#define EMIT(...)                                                              \
  if (sisi4s::Emitter::getRank() != 0) {                                       \
  } else sisi4s::Emitter::getEmitter()
#define EMIT_FLUSH(...)                                                        \
  if (sisi4s::Emitter::getRank() != 0) {                                       \
  } else sisi4s::Emitter::flush()

#endif
