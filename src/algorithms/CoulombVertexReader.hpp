#ifndef COULOMB_VERTEX_READER_DEFINED
#define COULOMB_VERTEX_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/Tensor.hpp>
#include <cstdint>
#include <fstream>

namespace sisi4s {
/**
 * \brief Reads the Coulomb vertex \f$\Gamma_{pG}^q\f$ and the occupied and
 * virtual orbital energies \f$\varepsilon_i, \varepsilon_a\f$ from binary
 * data file, and stores them in the CTF Tensors GammaGqr, epsi, epsa.
 */
DEFINE_ALGORITHM_HEADER(

    CoulombVertexReader,

    class Header {
    public:
      char magic[8];
      int32_t No, Nv, NG, NSpins, kPoints, reserved_;
      static char const *MAGIC;
    };
    class Chunk {
    public:
      char magic[8];
      int64_t size;
      static char const *REALS_MAGIC;
      static char const *IMAGS_MAGIC;
      static char const *REALSIA_MAGIC;
      static char const *IMAGSIA_MAGIC;
      static char const *EPSILONS_MAGIC;
    };

    void handleUnrestricted();
    void unrestrictVertex();
    void unrestrictEigenEnergies(const std::string &name););
} // namespace sisi4s

#endif
