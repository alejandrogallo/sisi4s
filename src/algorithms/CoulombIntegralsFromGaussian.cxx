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
#include <Eigen/Eigenvalues>
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
  CoulombIntegralsProvider
    ( const size_t Np_
    , const libint2::BasisSet& shells_
    , const Operator op_ = Operator::coulomb
    ): Np(Np_) , shells(shells_) , op(op_)
   {

   const size_t np(Cc4s::world->np), rank(Cc4s::world->rank);
   const auto toCtfIdx([&](const size_t i) { return i * Np*Np*Np; });

   const size_t chunks = shells.size() / np;
   if (chunks == 0) throw EXCEPTION("Number of shells is less than nprocs");
   firstShellIdx = rank * chunks;
   lastShellIdx = (np - 1) == rank
                ? shells.size() - 1
                : (rank+1) * chunks - 1
                ;
   const ShellInfo first(shells, firstShellIdx), second(shells, lastShellIdx);
   beginIndex = toCtfIdx(first.begin);
   endIndex = toCtfIdx(second.end) - 1;
   LOGGER(1) << "chunks       :" << chunks << std::endl;
   LOGGER(1) << "firstShellIdx:" << firstShellIdx << std::endl;
   LOGGER(1) << "lastShellIdx :" << lastShellIdx << std::endl;
   LOGGER(1) << "beginIndex   :" << beginIndex << std::endl;
   LOGGER(1) << "endIndex     :" << endIndex << std::endl;
#ifdef DEBUG
   std::cout << rank << "::firstShellIdx:" << firstShellIdx << std::endl;
   std::cout << rank << "::lastShellIdx :" << lastShellIdx << std::endl;
   std::cout << rank << "::beginIndex   :" << beginIndex << std::endl;
   std::cout << rank << "::endIndex     :" << endIndex << std::endl;
   std::cout << rank << "::DeltaIndex   :" << endIndex - beginIndex << std::endl;
#endif
  }

  void compute() {
    libint2::initialize();

    const size_t NpNpNpNp(Np*Np*Np*Np)
               , localSize(size())
               ;
    LOGGER(1) << "Np**4   :  " << NpNpNpNp << std::endl;
    LOGGER(1) << " | local:  " << localSize << std::endl;
    LOGGER(1) << "Allocating and computing Vklmn ("
              << sizeof(double) * NpNpNpNp / 1024 / 1024 / 1024
              << " GB)" << std::endl;
    Vklmn.resize(localSize, 0.0);

    libint2::Engine engine(op, shells.max_nprim(), shells.max_l(), 0);

    // store shell by shell calculation in this buffer
    const auto& vsrqp = engine.results();

    // the outside loops will loop over the shells.
    // This will create a block of Vpqrs, where pqrs are contracted
    // gaussian indices belonging to their respective shells.
    for (size_t _K(firstShellIdx); _K < lastShellIdx+1; ++_K) { // kappa
    for (size_t _L(0); _L < shells.size(); ++_L) { // lambda
    for (size_t _M(0); _M < shells.size(); ++_M) { // mu
    for (size_t _N(0); _N < shells.size(); ++_N) { // nu
      const ShellInfo K(shells, _K)
                    , L(shells, _L)
                    , M(shells, _M)
                    , N(shells, _N)
                    ;

      // compute integrals (K L , M N)
      engine.compute(shells[_K], shells[_L], shells[_M], shells[_N]);

        //for (size_t k(K.begin), Inmlk = 0; k < K.end; ++k) {
        for (size_t k(K.begin), Inmlk = 0; k < K.end; ++k) {
        for (size_t l(L.begin); l < L.end; ++l) {
        for (size_t m(M.begin); m < M.end; ++m) {
        for (size_t n(N.begin); n < N.end; ++n, ++Inmlk) {

          const size_t bigI( n
                           + m * Np
                           + l * Np*Np
                           + k * Np*Np*Np
                           );

          Vklmn[bigI - beginIndex] += vsrqp[0][Inmlk];

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

  double* data() {return Vklmn.data();}
  size_t size() const { return end() - begin() + 1; }
  size_t begin() const { return beginIndex; }
  size_t end() const { return endIndex; }

  private:
  const size_t Np;
  const libint2::BasisSet& shells;
  const Operator op;
  std::vector<double> Vklmn;
  size_t firstShellIdx, lastShellIdx, beginIndex, endIndex;
};


void CoulombIntegralsFromGaussian::run() {

  const std::vector<std::string> allArguments =
    { "xyzStructureFile"
    , "basisSet"
    , "kernel"
    , "chemistNotation"
    , "CoulombIntegrals"
    };
  checkArgumentsOrDie(allArguments);

  const std::string xyzStructureFile(getTextArgument("xyzStructureFile", ""))
                  , basisSet(getTextArgument("basisSet"))
                  , kernel(getTextArgument("kernel", "coulomb"))
                  ;
  const bool chemistNotation(getIntegerArgument("chemistNotation", 1) == 1);

  std::ifstream structureFileStream(xyzStructureFile.c_str());
  const auto atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();
  const libint2::BasisSet shells(basisSet, atoms);
  const int Np(shells.nbf());
  const Operator op([&](void){
    if      (kernel == "coulomb") return Operator::coulomb;
    else if (kernel == "delta")   return Operator::delta;
    else                          return Operator::coulomb;
  }());

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

  CoulombIntegralsProvider engine(Np, shells, op);

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
  {
    std::vector<int64_t> indices(engine.size());
    // important to start at engine.begin()
    std::iota(indices.begin(), indices.end(), engine.begin());
#ifdef DEBUG
    std::cout << Cc4s::world->rank
              << "::writing " << indices.size() << " values "
              << "into ctf tensor" << std::endl;
#endif
    LOGGER(1) << "writing " << indices.size() << " values "
              << "into ctf tensor" << std::endl;
    Vklmn->write(indices.size(), indices.data(), engine.data());
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
