#ifndef LINEAR_MIXER_DEFINED
#define LINEAR_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace sisi4s {
template <typename F>
class LinearMixer : public Mixer<F> {
public:
  MIXER_REGISTRAR_DECLARATION(LinearMixer);
  LinearMixer(Algorithm *algorithm);
  virtual ~LinearMixer();

  virtual void append(const PTR(FockVector<F>) &A, const PTR(FockVector<F>) &R);
  virtual PTR(const FockVector<F>) get();
  virtual PTR(const FockVector<F>) getResiduum();

  PTR(FockVector<F>) last;
  PTR(FockVector<F>) lastResiduum;

  double ratio;
};
} // namespace sisi4s

#endif
