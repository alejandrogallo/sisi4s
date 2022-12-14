#ifndef MPI_COMMUNICATOR_DEFINED
#define MPI_COMMUNICATOR_DEFINED

#include "mpi.h"
#include <math/Complex.hpp>
#include <math/Vector.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
// base template for type traits
template <typename F>
class MpiTypeTraits;

class MpiCommunicator {
public:
  MpiCommunicator(int rank_, int processes_, MPI_Comm comm_ = MPI_COMM_WORLD)
      : rank(rank_)
      , processes(processes_)
      , comm(comm_) {}
  MpiCommunicator(const CTF::World &world)
      : rank(world.rank)
      , processes(world.np)
      , comm(world.comm) {}
  ~MpiCommunicator() {}

  void barrier() { MPI_Barrier(comm); }

  template <typename F>
  void reduce(const F &src, F &dst, int rootRank = 0) {
    MPI_Reduce(&src,
               &dst,
               MpiTypeTraits<F>::elementCount(),
               MpiTypeTraits<F>::elementType(),
               MPI_SUM,
               rootRank,
               comm);
  }

  template <typename F>
  void allReduce(const F &src, F &dst) {
    MPI_Allreduce(&src,
                  &dst,
                  MpiTypeTraits<F>::elementCount(),
                  MpiTypeTraits<F>::elementType(),
                  MPI_SUM,
                  comm);
  }

  /**
   * \Brief Gathers the src vectors of all ranks together to the dst
   * vector at the given root rank, by default rank 0.
   **/
  template <typename F>
  void
  gather(const std::vector<F> &src, std::vector<F> &dst, int rootRank = 0) {
    if (rank == rootRank) {
      dst.resize(src.size() * processes);
    } else {
      dst.resize(0);
    }
    MPI_Gather(src.data(),
               src.size() * MpiTypeTraits<F>::elementCount(),
               MpiTypeTraits<F>::elementType(),
               dst.data(),
               src.size() * MpiTypeTraits<F>::elementCount(),
               MpiTypeTraits<F>::elementType(),
               rootRank,
               comm);
  }

  int getRank() const { return rank; }
  int getProcesses() const { return processes; }

protected:
  int rank, processes;
  MPI_Comm comm;
};

template <>
class MpiTypeTraits<int> {
public:
  static MPI_Datatype elementType() { return MPI_INT; }
  static int elementCount() { return 1; }
};

template <>
class MpiTypeTraits<int64_t> {
public:
  static MPI_Datatype elementType() { return MPI_INTEGER8; }
  static int elementCount() { return 1; }
};

template <>
class MpiTypeTraits<uint64_t> {
public:
  static MPI_Datatype elementType() { return MPI_INTEGER8; }
  static int elementCount() { return 1; }
};

template <>
class MpiTypeTraits<double> {
public:
  static MPI_Datatype elementType() { return MPI_REAL8; }
  static int elementCount() { return 1; }
};

template <>
class MpiTypeTraits<complex> {
public:
  static MPI_Datatype elementType() { return MPI_DOUBLE_COMPLEX; }
  static int elementCount() { return 1; }
};

template <typename F, int D>
class MpiTypeTraits<Vector<F, D>> {
public:
  static MPI_Datatype elementType() { return MpiTypeTraits<F>::elementType(); }
  static int elementCount() { return D; }
};
} // namespace sisi4s

#endif
