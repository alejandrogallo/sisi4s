#include <string>
#include <vector>
#include <algorithm>
#include <libint2.hpp>
#include <algorithms/CoulombIntegralsFromGaussian.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#include <set>
#include <map>
#include <util/Emitter.hpp>
#define LOGGER(_l) LOG(_l, "CoulombIntegralsFromGaussian")

typedef libint2::Operator Operator;

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromGaussian);

// struct for storing information about the shell ends in the for loops
// calculating the integrals
struct ShellInfo {
  // size of the shell, global begin and global end
  // l should be the angular momentum
  size_t size, begin, end, l;
  // constructor from a BasisSet and a shell index i
  inline ShellInfo(const libint2::BasisSet &shells, const size_t i) {
    size = shells[i].size();
    begin = shells.shell2bf()[i];
    end = size + begin;
    l = shells[i].contr[0].l;
  }
};

struct CoulombIntegralsProvider {
  enum Distribution { SIMPLE, ROUND_ROBIN };
  CoulombIntegralsProvider( const libint2::BasisSet& shells_
                          , const Operator op_
                          , const Distribution d
                          ) : Np(shells_.nbf()), shells(shells_), op(op_)
  {
    // set mpiShells
    if (d == SIMPLE) doSimpleDistribution();
    else             doRoundRobinDistribution();

    const size_t np(Cc4s::world->np);
    // this we need to be able to tell every np how many times
    // it should perform an mpi operation in the case it is necessary
    maxSize = (shells.size() < np)
            ? 1
            : shells.size() / np + shells.size() % np
            ;

    LOGGER(1) << "maxSize: " << maxSize << std::endl;

    { std::stringstream s; // print the mpiShells
      for (const auto& i: mpiShells) s << i << ",";
      LOGGER(1) << "shells: {" << s.str() << "}" << std::endl;
    }

    const auto toCtfIdx([&](const size_t i) { return i * Np*Np*Np; });
    // to store information about the shells
    std::vector<ShellInfo> _shells;

    // fill in the shells information
    std::transform( mpiShells.begin()
                  , mpiShells.end()
                  , std::back_inserter(_shells)
                  , [&](const size_t &i) { return ShellInfo(shells, i);}
                  );

    // get the ctf indices
    for (auto &s: _shells) {
      const size_t bx(toCtfIdx(s.begin)), ex(toCtfIdx(s.end));
      // LOGGER(1) << "sb: "<< s.begin << " se: "<< s.end << std::endl;
      // LOGGER(1) << "bx: "<< bx << " ex: "<< ex << std::endl;
      std::vector<int64_t> idxctf(ex - bx);
      std::iota(idxctf.begin(), idxctf.end(), toCtfIdx(s.begin));
      ctfIndices.push_back(idxctf);
    }

    LOGGER(1) << "#shells        :" << mpiShells.size()  << std::endl;
    LOGGER(1) << "#indices       :" << ctfIndices.size() << std::endl;

  }

  void doRoundRobinDistribution() {
    const size_t np(Cc4s::world->np), rank(Cc4s::world->rank);

    for (size_t s(0); s < shells.size(); s++)
      if (s % np == rank)
        mpiShells.push_back(s);

  }

  void doSimpleDistribution() {
    const size_t np   = Cc4s::world->np
               , rank = Cc4s::world->rank
               , np_eff = std::min(np, shells.size())
               , chunks = shells.size() / np_eff
               , firstShellIdx = rank * chunks
               , lastShellIdx  =
                   [&] {
                     if ((np_eff - 1) == rank) return shells.size();
                     else if (rank >= np_eff)  return firstShellIdx;
                     else                      return (rank+1) * chunks;
                   }()
               ;
    mpiShells.resize(lastShellIdx - firstShellIdx);
    std::iota(mpiShells.begin(), mpiShells.end(), firstShellIdx);

    LOGGER(1) << "np effective " << np_eff << std::endl;
    LOGGER(1) << "chunks       :" << chunks << std::endl;
    LOGGER(1) << "firstShellIdx:" << firstShellIdx << std::endl;
    LOGGER(1) << "lastShellIdx :" << lastShellIdx << std::endl;
#ifdef DEBUGG
    std::cout << rank << "::firstShellIdx:" << firstShellIdx << std::endl;
    std::cout << rank << "::lastShellIdx :" << lastShellIdx << std::endl;
#endif
  }

  void compute() {
    libint2::initialize();

    const size_t NpNpNpNp(Np*Np*Np*Np);
    LOGGER(1) << "Np**4   :  " << NpNpNpNp << std::endl;
    LOGGER(1) << "Computing Vklmn ("
              << sizeof(double) * NpNpNpNp / 1024 / 1024 / 1024
              << " GB)" << std::endl;

    libint2::Engine engine( op
                          , shells.max_nprim()
                          , shells.max_l()
                          , 0
                          );

    // store shell by shell calculation in this buffer
    const auto& vsrqp = engine.results();

    // the outside loops will loop over the shells.
    // This will create a block of Vpqrs, where pqrs are contracted
    // gaussian indices belonging to their respective shells.
    int mpiShellIdx(-1);
    for (const size_t _K: mpiShells) {
      ++mpiShellIdx;
      LOGGER(1) << "shell: " << mpiShells[mpiShellIdx]
                << "    #: " << ctfIndices[mpiShellIdx].size()
                << std::endl;
      Vklmn.push_back(std::vector<double>(ctfIndices[mpiShellIdx].size(), 0));
    for (size_t _L(0); _L < shells.size(); ++_L         ) { // lambda
    for (size_t _M(0); _M < shells.size(); ++_M         ) { // mu
    for (size_t _N(0); _N < shells.size(); ++_N         ) { // nu
      const ShellInfo K(shells, _K)
                    , L(shells, _L)
                    , M(shells, _M)
                    , N(shells, _N)
                    ;

      // compute integrals (K L , M N)
      engine.compute(shells[_K], shells[_L], shells[_M], shells[_N]);

      if (vsrqp[0] == nullptr) continue;

      for (size_t k(K.begin), Inmlk = 0; k < K.end; ++k         ) {
      for (size_t l(L.begin)           ; l < L.end; ++l         ) {
      for (size_t m(M.begin)           ; m < M.end; ++m         ) {
      for (size_t n(N.begin)           ; n < N.end; ++n, ++Inmlk) {

        const size_t bigI( n
                         + m * Np
                         + l * Np*Np
                         + k * Np*Np*Np
                         )
                   , idx(bigI - ctfIndices[mpiShellIdx].front())
                   ;

        Vklmn.back()[idx] += vsrqp[0][Inmlk];

      } // n
      } // m
      } // l
      } // k

    } // N
    } // M
    } // L
    } // K

  libint2::finalize();

  }

  std::vector<std::vector<double>>     &data() { return Vklmn; }
  std::vector<std::vector<int64_t>> &indices() { return ctfIndices; }
  // mpi max size of Vklmn and ctfIndices
  // to be able to write Kosher with regard to mpi
  size_t maxSize;

  private:
  const size_t Np;
  const libint2::BasisSet& shells;
  const Operator op;
  std::vector<size_t> mpiShells;
  std::vector<std::vector<double>> Vklmn;
  std::vector<std::vector<int64_t>> ctfIndices;
};


void CoulombIntegralsFromGaussian::run() {

  checkArgumentsOrDie( { "xyzStructureFile"
                       , "basisSet"
                       , "kernel"
                       , "chemistNotation"
                       , "shellDistribution"
                       , "CoulombIntegrals"
                       } );

  const std::string xyzStructureFile(getTextArgument("xyzStructureFile", ""))
                  , basisSet(getTextArgument("basisSet"))
                  , kernel(getTextArgument("kernel", "coulomb"))
                  ;
  const bool chemistNotation(getIntegerArgument("chemistNotation", 1) == 1);
  const CoulombIntegralsProvider::Distribution
    mpiDistribution = [&]{
      auto mode(this->getTextArgument("shellDistribution", "simple"));
      CoulombIntegralsProvider::Distribution r;
      if      (mode == "simple")     r = CoulombIntegralsProvider::SIMPLE;
      else if (mode == "roundRobin") r = CoulombIntegralsProvider::ROUND_ROBIN;
      else                           throw "Incorrect mode given: " + mode;
      LOGGER(1) << "mpi mode: " << mode << " (" << r << ")" << std::endl;
      return r;
    }();

  std::ifstream structureFileStream(xyzStructureFile);
  if (!structureFileStream.good()) throw "Bad file: " + xyzStructureFile;
  const auto atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();
  const libint2::BasisSet shells(basisSet, atoms);
  const int Np(shells.nbf());
  const Operator
    op = [kernel]{ if      (kernel == "coulomb") return Operator::coulomb;
                   else if (kernel == "delta")   return Operator::delta;
                   else                          throw "Operator not valid";
                 }();

  LOGGER(1) << "libint          : " << LIBINT_VERSION << std::endl;
  LOGGER(1) << "libint max_am   : " << LIBINT_MAX_AM << std::endl;
  LOGGER(1) << "cgshell_ordering: " << LIBINT_CGSHELL_ORDERING << std::endl;
  LOGGER(1) << "std ordering is : " << LIBINT_CGSHELL_ORDERING_STANDARD
            << std::endl;
  LOGGER(1) << "Np: " << Np << std::endl;
  LOGGER(1) << "kernel: " << kernel << std::endl;
  LOGGER(1) << "structure: " << xyzStructureFile << std::endl;
  LOGGER(1) << "basisSet: " << basisSet << std::endl;
  LOGGER(1) << "#shells: " << shells.size() << std::endl;

  CoulombIntegralsProvider engine(shells, op, mpiDistribution);

  const std::vector<int> lens(4, Np);
  const std::vector<int> syms(4, NS);
  auto Vklmn(new CTF::Tensor<double>(4, lens.data(), syms.data(),
                                      *Cc4s::world, "V"));

  EMIT() << YAML::Key << "Np" << YAML::Value << Np
         << YAML::Key << "xyz" << YAML::Value << xyzStructureFile
         << YAML::Key << "basis-set" << YAML::Value << basisSet
         << YAML::Key << "kernel" << YAML::Value << kernel
         << YAML::Key << "chemist-notation" << YAML::Value << chemistNotation
         ;

  engine.compute();
  LOGGER(1) << "computation done" << std::endl;

  // loop over engine.maxSize to ensure the same number of mpi calls
  // on every rank
  for (size_t i(0); i < engine.maxSize; i++) {
    int64_t* indices;
    size_t N;
    double* data;
    if (i < engine.indices().size()) {
      N = engine.indices()[i].size();
      indices = engine.indices()[i].data();
      data = engine.data()[i].data();
    } else {
      N = 0;
      indices = nullptr;
      data = nullptr;
    }
#ifdef DEBUG
    std::cout << Cc4s::world->rank << "::"
              << "writing " << N << " values "
              << "into ctf tensor" << std::endl;
#endif
    LOGGER(1) << "writing " << N << " values "
              << "into ctf tensor" << std::endl;
    Vklmn->write(N, indices, data);
  }

  if (chemistNotation) {

    LOGGER(1) << "WARNING: this integral is in chemist notation,"
                 " Vklmn = (kl|mn) = <km|ln>" << std::endl;
    allocatedTensorArgument<double>("CoulombIntegrals", Vklmn);

  } else {

    LOGGER(1) << "Allocating phyisics notation integrals" << std::endl;
    auto newV(new CTF::Tensor<double>(4, Vklmn->lens, Vklmn->sym, *Cc4s::world));
    LOGGER(1) << "Converting chemist integral into physics" << std::endl;
    LOGGER(1) << "| PhysV[pqrs] = ChemV[prqs]; " << std::endl;

    (*newV)["pqrs"] = (*Vklmn)["prqs"];
    LOGGER(1) << "NOTE: this integral is in physicist notation" << std::endl;
    allocatedTensorArgument<double>("CoulombIntegrals", newV);
    delete Vklmn;

  }

}
