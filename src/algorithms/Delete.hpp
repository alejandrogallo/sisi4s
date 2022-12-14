#ifndef DELETE_DEFINED
#define DELETE_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class Delete : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(Delete);
  Delete(std::vector<Argument> const &argumentList);
  virtual ~Delete();
  /**
   * \brief Frees all associated resources of the given tensor
   */
  virtual void run();
  /**
   * \brief Frees all associated dry resources of the given tensor
   */
  virtual void dryRun();
};
} // namespace sisi4s

#endif
