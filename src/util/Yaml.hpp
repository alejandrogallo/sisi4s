#include <yaml-cpp/yaml.h>

#define YAML_ASSERT_KEY(node, key)                                             \
  do {                                                                         \
    if (!node[key].IsDefined()) {                                              \
      std::stringstream s;                                                     \
      const auto mark = node.Mark();                                           \
      s << "Node in line:col " << mark.line << ":" << mark.column << "\n\n"    \
        << node << "\n\n"                                                      \
        << "should have a key " << key << "\n\n";                              \
      throw s.str();                                                           \
    }                                                                          \
  } while (0)
