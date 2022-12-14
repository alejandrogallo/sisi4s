#ifndef CTF_MACHINE_TENSOR_DEFINED
#define CTF_MACHINE_TENSOR_DEFINED

#include <tcc/MachineTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <string>
#include <memory>

namespace sisi4s {
template <typename F = double>
class CtfMachineTensorFactory;

template <typename F = double>
class CtfMachineTensor : public tcc::MachineTensor<F> {
protected:
  class ProtectedToken {};

public:
  // required by templates to infer corresponding Factory type
  typedef CtfMachineTensorFactory<F> Factory;
  typedef Tensor<F> Tensor;

  // constructors called by factory
  CtfMachineTensor(const std::vector<int> &lens,
                   const std::string &name,
                   CTF::World *world,
                   const ProtectedToken &)
      : tensor(static_cast<int>(lens.size()),
               lens.data(),
               std::vector<int>(0, lens.size()).data(),
               *world,
               name.c_str()) {}

  // copy constructor from CTF tensor, for compatibility
  CtfMachineTensor(const Tensor &T, const ProtectedToken &)
      : tensor(T) {}

  static std::shared_ptr<CtfMachineTensor<F>> create(const Tensor &T) {
    return std::make_shared<CtfMachineTensor<F>>(T, ProtectedToken());
  }

  virtual ~CtfMachineTensor() {}

  // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
  virtual void move(F alpha,
                    const std::shared_ptr<tcc::MachineTensor<F>> &A,
                    const std::string &aIndices,
                    F beta,
                    const std::string &bIndices) {
    std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A));
    if (!ctfA) {
      throw new EXCEPTION("Passed machine tensor of wrong implementation.");
    }
    LOG(2, "TCC") << "move " << getName() << "[" << bIndices
                  << "] <<= " << alpha << " * " << ctfA->getName() << "["
                  << aIndices << "] + " << beta << " * " << getName() << "["
                  << bIndices << "]" << std::endl;
    tensor.sum(alpha, ctfA->tensor, aIndices.c_str(), beta, bIndices.c_str());
  }

  // this[bIndices] = alpha * f(A[aIndices]) + beta*this[bIndices]
  void move(F alpha,
            const std::shared_ptr<tcc::MachineTensor<F>> &A,
            const std::string &aIndices,
            F beta,
            const std::string &bIndices,
            const std::function<F(const F)> &f) {
    std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A));
    if (!ctfA) {
      throw new EXCEPTION("Passed machine tensor of wrong implementation.");
    }
    LOG(2, "TCC") << "move " << getName() << "[" << bIndices
                  << "] <<= " << alpha << " * " << ctfA->getName() << "["
                  << aIndices << "] + " << beta << " * " << getName() << "["
                  << bIndices << "]" << std::endl;
    tensor.sum(alpha,
               ctfA->tensor,
               aIndices.c_str(),
               beta,
               bIndices.c_str(),
               CTF::Univar_Function<F>(f));
  }

  // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
  void contract(F alpha,
                const std::shared_ptr<tcc::MachineTensor<F>> &A,
                const std::string &aIndices,
                const std::shared_ptr<tcc::MachineTensor<F>> &B,
                const std::string &bIndices,
                F beta,
                const std::string &cIndices) {
    std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A));
    std::shared_ptr<CtfMachineTensor<F>> ctfB(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(B));
    if (!ctfA || !ctfB) {
      throw new EXCEPTION("Passed machine tensor of wrong implementation.");
    }
    LOG(2, "TCC") << "contract " << getName() << "[" << cIndices << "] <<= g("
                  << alpha << " * " << ctfA->getName() << "[" << aIndices
                  << "], " << ctfB->getName() << "[" << bIndices << "]) + "
                  << beta << " * " << getName() << "[" << cIndices << "]"
                  << std::endl;
    tensor.contract(alpha,
                    ctfA->tensor,
                    aIndices.c_str(),
                    ctfB->tensor,
                    bIndices.c_str(),
                    beta,
                    cIndices.c_str());
  }

  // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
  void contract(F alpha,
                const std::shared_ptr<tcc::MachineTensor<F>> &A,
                const std::string &aIndices,
                const std::shared_ptr<tcc::MachineTensor<F>> &B,
                const std::string &bIndices,
                F beta,
                const std::string &cIndices,
                const std::function<F(const F, const F)> &g) {
    std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A));
    std::shared_ptr<CtfMachineTensor<F>> ctfB(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(B));
    if (!ctfA || !ctfB) {
      throw new EXCEPTION("Passed machine tensor of wrong implementation.");
    }
    LOG(2, "TCC") << "contract " << getName() << "[" << cIndices << "] <<= g("
                  << alpha << " * " << ctfA->getName() << "[" << aIndices
                  << "], " << ctfB->getName() << "[" << bIndices << "]) + "
                  << beta << " * " << getName() << "[" << cIndices << "]"
                  << std::endl;
    tensor.contract(alpha,
                    ctfA->tensor,
                    aIndices.c_str(),
                    ctfB->tensor,
                    bIndices.c_str(),
                    beta,
                    cIndices.c_str(),
                    CTF::Bivar_Function<F>(g));
  }

  // TODO: interfaces to be defined: slice, permute, transform

  virtual std::vector<int> getLens() const {
    return std::vector<int>(tensor.lens, tensor.lens + tensor.order);
  }

  virtual std::string getName() const { return std::string(tensor.get_name()); }

  /**
   * \brief The adapted CTF tensor
   **/
  Tensor tensor;

  friend class CtfMachineTensorFactory<F>;
};

template <typename F>
class CtfMachineTensorFactory : public tcc::MachineTensorFactory<F> {
protected:
  class ProtectedToken {};

public:
  CtfMachineTensorFactory(CTF::World *world_, const ProtectedToken &)
      : world(world_) {}

  virtual ~CtfMachineTensorFactory() {}

  virtual std::shared_ptr<tcc::MachineTensor<F>>
  createTensor(const std::vector<int> &lens, const std::string &name) {
    return std::shared_ptr<typename tcc::MachineTensor<F>>(
        std::make_shared<CtfMachineTensor<F>>(
            lens,
            name,
            world,
            typename CtfMachineTensor<F>::ProtectedToken()));
  }

  static std::shared_ptr<CtfMachineTensorFactory<F>>
  create(CTF::World *world = Sisi4s::world) {
    return std::make_shared<CtfMachineTensorFactory<F>>(world,
                                                        ProtectedToken());
  }

protected:
  CTF::World *world;
};
} // namespace sisi4s

#endif
