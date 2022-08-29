#include <string>
#include <vector>
#include <algorithm>
#include <libint2.hpp>
#include <algorithms/CoulombIntegralsFromGaussian.hpp>
#include <util/CTF.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <util/CTF.hpp>
#include <numeric>
#include <set>
#include <map>
#include <util/Emitter.hpp>
#include <math/MathFunctions.hpp>

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
                          , const size_t maxElementsWrite_
                          )
                          : Np(shells_.nbf())
                          , shells(shells_)
                          , op(op_)
                          , maxElementsWrite(maxElementsWrite_)
                          , maxShellNbf(2*shells.max_l() + 1)
                          , maxShellElements(maxShellNbf * Np*Np*Np)
                          , writesPerShell( maxShellElements
                                          / maxElementsWrite + 1
                                          )
  {

    const auto toCtfIdx([&](const size_t i) { return i * Np*Np*Np; });

    // set mpiShells
    if (d == SIMPLE) doSimpleDistribution();
    else             doRoundRobinDistribution();

    // How many writes will we make per shells is controlled by the
    // ratio of maxShellElements / maxElementsWrite
    LOGGER(1) << "writes per shell   : " << writesPerShell   << std::endl;
    LOGGER(1) << "max shell l        : " << shells.max_l()   << std::endl;
    LOGGER(1) << "max shell nbf      : " << maxShellNbf      << std::endl;
    LOGGER(1) << "max shell elements : " << maxShellElements << std::endl;
    LOGGER(1) << "max elem write     : " << maxElementsWrite << std::endl;

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

  void computeHalfContract(size_t No, CTF::Tensor<double> *coeff) {
    libint2::initialize();

    size_t orbs(coeff->lens[1]);
    LOGGER(1) << "dimensions: " << coeff->lens[0] << " " << coeff->lens[1] << std::endl;
    if ( coeff->lens[0] != Np ) throw "dimension problems!";

    std::vector<double> C(orbs*Np);
    coeff->read_all(C.data());
    const size_t NoNoNpNp(No*No*Np*Np);
    LOGGER(1) << "No**2*Np**2   :  " << NoNoNpNp << std::endl;
    LOGGER(1) << "Computing half transformed Vkimj ("
              << (double) sizeof(double) * NoNoNpNp / 1024 / 1024 / 1024
              << " GB)" << std::endl;
    LOGGER(1) << "Largest buffer size for indiviual mpi rank: " <<
                 (double) sizeof(double) * maxShellElements /1024/1024/1024 << " GB\n";
    libint2::Engine engine( op
                          , shells.max_nprim()
                          , shells.max_l()
                          , 0
                          );

    // store shell by shell calculation in this buffer
    const auto& vsrqp = engine.results();
    // store the values Vklmn values in pppp
    // half-transform pppp -> pppo -> ppoo
    std::vector<double> pppp(maxShellElements);
    std::vector<double> pppo(maxShellNbf*Np*Np*No);
    std::vector<double> ppoo(maxShellNbf*Np*No*No);
    offset.resize(ctfIndices.size(),-1);
    // the outside loops will loop over the shells.
    // This will create a block of Vpqrs, where pqrs are contracted
    // gaussian indices belonging to their respective shells.
    int mpiShellIdx(-1);
    for (const size_t _K: mpiShells) {
      ++mpiShellIdx;
      LOGGER(1) << "shell: " << mpiShells[mpiShellIdx]
                << "    #: " << ctfIndices[mpiShellIdx].size()
                << std::endl;
      const ShellInfo K(shells, _K);
      offset[mpiShellIdx] = K.begin*Np*No*No;
//      std::cout << "offset " << K.begin*Np*No*No << std::endl;
      Vklmn.push_back(std::vector<double>(K.size*Np*No*No, 0.));
      std::fill(pppp.begin(), pppp.end(), 0.);
      std::fill(pppo.begin(), pppo.end(), 0.);
      std::fill(ppoo.begin(), ppoo.end(), 0.);
//      Vklmn.push_back(std::vector<double>(ctfIndices[mpiShellIdx].size(), 0));
    for (size_t _L(0); _L < shells.size(); ++_L         ) { // lambda
    for (size_t _M(0); _M < shells.size(); ++_M         ) { // mu
    for (size_t _N(0); _N < shells.size(); ++_N         ) { // nu
      const ShellInfo L(shells, _L)
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

        pppp[idx] += vsrqp[0][Inmlk];

      } // n
      } // m
      } // l
      } // k

    } // N
    } // M
    } // L
      // Now we finished writing Vklmn, i.e. the AO Coulomb Matrix for the shell K
      // We can now contract Vklmn.back() with the Orbital coefficient C_li C_nj


      // contract:    V(klmi) = V(klmn) * C(ni);
      for (size_t k(0); k < K.size; k++)
      for (size_t l(0); l < Np;     l++)
      for (size_t m(0); m < Np;     m++){
        size_t offi(k*Np*Np*Np + l*Np*Np + m*Np);
        size_t offo(k*Np*Np*No + l*Np*No + m*No);
      for (size_t n(0); n < Np;     n++){
      for (size_t i(0); i < No;     i++){
        pppo[offo + i] += pppp[offi + n] * C[i*Np + n];
      } // i
      } // n
      } // m

      // rearrange pppo: V(klmi) -> V(kmil)
      // this means: 1.) particle index l getting the fastest index
      //             2.) chemists integral changes to physics notation
      // we use pppp as a scratch array
      for (size_t s(0); s < K.size*Np*No*Np; s++)  pppp[s] = pppo[s];
      std::fill(pppo.begin(), pppo.end(), 0.);

      for (size_t k(0); k < K.size; k++)
      for (size_t m(0); m < Np;     m++)
      for (size_t i(0); i < No;     i++)
      for (size_t l(0); l < Np;     l++)
        pppo[k*Np*No*Np + m*Np*No + i*Np + l] = pppp[k*Np*No*Np + l*Np*No + m*No + i];



      // contract: V(kmij) = V(kmil) * C(lj)
      for (size_t k(0); k < K.size; k++)
      for (size_t m(0); m < Np;     m++)
      for (size_t i(0); i < No;     i++){
        size_t offi(k*Np*Np*No + m*Np*No + i*Np);
        size_t offo(k*Np*No*No + m*No*No + i*No);
      for (size_t l(0); l < Np;     l++){
      for (size_t j(0); j < No;     j++){
        ppoo[offo + j] += pppo[offi + l] * C[j*Np + l];
      } // j
      } // l
      } // i

      for (size_t s(0); s < K.size*Np*No*No; s++) Vklmn.back()[s] = ppoo[s];



    } // K

  libint2::finalize();

  }


  struct IndexRange {
    size_t begin, end, size;
    IndexRange(): begin(0), end(0), size(0) {};
  };
  std::vector<IndexRange> getWriteRangesForShell(size_t i) {

    std::vector<IndexRange> r(writesPerShell);

    if (i >= ctfIndices.size()) return r;
    const auto& indices(ctfIndices[i]);

    const size_t chunks = indices.size() / maxElementsWrite
               , rest = indices.size() % maxElementsWrite
               ;

    size_t j(0);
    for (j = 0; j < chunks; j++) {
      r[j].begin  = maxElementsWrite * j;
      r[j].end = maxElementsWrite * (j+1) - 1;
      r[j].size = r[j].end - r[j].begin + 1;
    }

    if (rest != 0) {
      r[j].begin  = maxElementsWrite * j;
      r[j].end = r[j].begin + rest - 1;
      r[j].size = r[j].end - r[j].begin + 1;
    }

    return std::move(r);

  }

  void write(CTF::Tensor<double> *t) {

    // loop over engine.maxSize to ensure the same number of mpi calls
    // on every rank
    for (size_t i(0); i < maxSize; i++) {
      int64_t* indices;
      size_t N;
      double* data;
      if (i < ctfIndices.size()) {
        N = ctfIndices[i].size();
        indices = ctfIndices[i].data();
        data = Vklmn[i].data();
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

      // ip is an index pair
      for (const auto &ip: getWriteRangesForShell(i)) {
        #ifdef DEBUG
        std::cout << Cc4s::world->rank << "::" << "\t|"
                  << "writing " << ip.size << " values "
                  << "into ctf tensor" << std::endl;
        #endif
        LOGGER(1) << "writing " << ip.size << " values "
                  << "into ctf tensor" << std::endl;
        t->write( ip.size
                , N ? indices + ip.begin : nullptr
                , N ? data    + ip.begin : nullptr
                );
      }

    }
  }

  void readHalfContract(CTF::Tensor<double> *t) {
    // loop over engine.maxSize to ensure the same number of mpi calls
    // on every rank
    // i.e. only one rank provide data for ctf->write.

    size_t off(0); int i(0);
    while (1) {
      size_t N;
      std::vector<int64_t> indices;
      double *data;
      size_t old(off);
      if ( (Vklmn.size() > 0) && (offset[i] == off)  ) {
        N = Vklmn[i].size();
        data = Vklmn[i].data();
        indices.resize(N);
        std::iota(indices.begin(), indices.end(), off);
        off += N;
        i++;
//        std::cout << "I am rank " << Cc4s::world->rank << " having offset " << off << std::endl;
      }
      else{
        N = 0;
      }
//      if (N>0) std::cout << "reading " << N << "elements\n";
      MPI_Allreduce(&off, &off, 1, MPI_INT64_T, MPI_MAX, Cc4s::world->comm);
      if ( off == old) break;
      t->write(N, N ? indices.data() : nullptr, N ? data : nullptr);
//      std::cout << "I am rank " << Cc4s::world->rank << " writing " << N
//                << " elements with offset " << off << std::endl;
    }



  }

  private:

    // mpi max size of Vklmn and ctfIndices
    // to be able to write Kosher with regard to mpi
    size_t maxSize;
    const size_t Np;
    const libint2::BasisSet& shells;
    const Operator op;
    std::vector<size_t> mpiShells;
    std::vector<size_t> offset;
    std::vector<std::vector<double>> Vklmn;
    std::vector<std::vector<int64_t>> ctfIndices;
    // this is 2**26, it is roughly 500 MB of double precission fp numbers
    const size_t maxElementsWrite;
    // maximum number of primitives in the basis set
    const size_t maxShellNbf;
    const size_t maxShellElements;
    const size_t writesPerShell;

};


void CoulombIntegralsFromGaussian::run() {

  checkArgumentsOrDie( { "xyzStructureFile"
                       , "basisSet"
                       , "kernel"
                       , "chemistNotation"
                       , "maxElementsWrite"
                       , "shellDistribution"
                       , "CoulombIntegrals"
                       , "HHHHCoulombIntegrals"
                       , "PPHHCoulombIntegrals"
                       , "OrbitalCoefficients"
                       , "HoleEigenEnergies"
                       , "Spins"
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
  const libint2::BasisSet shells(basisSet, atoms, true);
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

  EMIT() << YAML::Key << "Np" << YAML::Value << Np
         << YAML::Key << "xyz" << YAML::Value << xyzStructureFile
         << YAML::Key << "basis-set" << YAML::Value << basisSet
         << YAML::Key << "kernel" << YAML::Value << kernel
         << YAML::Key << "chemist-notation" << YAML::Value << chemistNotation
         ;


  CoulombIntegralsProvider engine( shells
                                 , op
                                 , mpiDistribution
                                 , getIntegerArgument( "maxElementsWrite"
                                                     , 67108864
                                                     )
                                 );

  if (isArgumentGiven("CoulombIntegrals")){


    const std::vector<int> lens(4, Np), syms(4, NS);
    auto Vklmn(new CTF::Tensor<double>( 4
                                      , lens.data()
                                      , syms.data()
                                      , *Cc4s::world, "V"
                                      ));

    LOGGER(1) << "computation begin" << std::endl;
    engine.compute();
    LOGGER(1) << "computation done" << std::endl;
    engine.write(Vklmn);
    LOGGER(1) << "writing done" << std::endl;

    if (chemistNotation) {

      LOGGER(1) << "WARNING: this integral is in chemist notation,"
                   " Vklmn = (kl|mn) = <km|ln>" << std::endl;
      allocatedTensorArgument<double>("CoulombIntegrals", Vklmn);

    } else {

      LOGGER(1) << "Allocating phyisics notation integrals" << std::endl;
      auto newV(new CTF::Tensor<double>( 4
                                       , Vklmn->lens
                                       , Vklmn->sym
                                       , *Cc4s::world
                                       ));
      LOGGER(1) << "Converting chemist integral into physics" << std::endl;
      LOGGER(1) << "| PhysV[pqrs] = ChemV[prqs]; " << std::endl;

      (*newV)["pqrs"] = (*Vklmn)["prqs"];
      LOGGER(1) << "NOTE: this integral is in physicist notation" << std::endl;
      allocatedTensorArgument<double>("CoulombIntegrals", newV);
      delete Vklmn;

    }
  }
  else if (  isArgumentGiven("OrbitalCoefficients")
          && isArgumentGiven("PPHHCoulombIntegrals")
          && isArgumentGiven("HoleEigenEnergies")) {

    auto epsi(getTensorArgument("HoleEigenEnergies"));
    const size_t No(epsi->lens[0]);

    auto C(getTensorArgument("OrbitalCoefficients"));
    const size_t Nv(C->lens[1] - No);
    std::vector<int> lens(4, Np), syms(4, NS);
    lens[0] = No; lens[1] = No;
    auto Vijkm(new CTF::Tensor<double>( 4, lens.data(), syms.data(), *Cc4s::world, "V"));

    LOGGER(1) << "computation begin" << std::endl;
    engine.computeHalfContract(No, C);
    LOGGER(1) << "computation done" << std::endl;
    engine.readHalfContract(Vijkm);
    LOGGER(1) << "ctfread done" << std::endl;




    lens[2] = No; lens[3] = No;
    auto Vhhhh(new CTF::Tensor<double>( 4, lens.data(), syms.data(), *Cc4s::world, "Vhhhh"));
    int sliceStart[] = {0, 0};
    int sliceEnd[] = {Np, No};
    auto Cocc(C->slice(sliceStart, sliceEnd));
    (*Vhhhh)["klij"] = (*Vijkm)["klpq"] * Cocc["qi"] * Cocc["pj"];
    if ( isArgumentGiven("Spins") ){
      LOGGER(1) << "unrestricted case: Vhhhh\n";
      auto S(getTensorArgument("Spins"));  
      int sS[] = {0}; int sE[] = {No};
      auto Socc(S->slice(sS, sE));
      std::vector<int> l(2);
      l[0] = No; l[1] = No;
      auto Sm(new CTF::Tensor<double>(2, l.data(), syms.data(), *Cc4s::world, "Sm"));
      (*Sm)["ij"] = Socc["i"]* Socc["j"];

      CTF::Transform<double>(
        std::function<void(double &)>(
           [](double &s) { s = (s + 0.25) * 2.0; }
        )
      )(
        (*Sm)["pq"]
      );
      auto Smap(new CTF::Tensor<double>(4, lens.data(), syms.data(), *Cc4s::world, "Smap"));
      (*Smap)["ijkl"] = (*Sm)["ik"] * (*Sm)["jl"];
      CTF::Bivar_Function<> fMultiply(&multiply<double>);
      Vhhhh->contract(
        1.0, *Vhhhh,"pqrs", *Smap,"pqrs", 0.0,"pqrs", fMultiply
      );
      delete Smap;
      delete Sm;
    }


    sliceStart[1] = No; sliceEnd[1] = C->lens[1];
    auto Cvirt(C->slice(sliceStart, sliceEnd));
    lens[0] = Nv; lens[1] = Nv;
    auto Vpphh(new CTF::Tensor<double>( 4, lens.data(), syms.data(), *Cc4s::world, "Vpphh"));
    (*Vpphh)["abij"] = Cvirt["qa"] * (*Vijkm)["ijpq"] * Cvirt["pb"];
    if ( isArgumentGiven("Spins") ){
      LOGGER(1) << "unrestricted case: Vpphh\n";
      auto S(getTensorArgument("Spins"));  
      int sS[] = {0}; int sE[] = {No};
      auto Socc(S->slice(sS, sE));
      sS[0] = No; sE[0] = C->lens[1];
      auto Svirt(S->slice(sS, sE));
      std::vector<int> l(2);
      l[0] = Nv; l[1] = No;
      auto Sm(new CTF::Tensor<double>(2, l.data(), syms.data(), *Cc4s::world, "Sm"));
      (*Sm)["ai"] = Svirt["a"]* Socc["i"];

      CTF::Transform<double>(
        std::function<void(double &)>(
           [](double &s) { s = (s + 0.25) * 2.0; }
        )
      )(
        (*Sm)["pq"]
      );
      auto Smap(new CTF::Tensor<double>(4, lens.data(), syms.data(), *Cc4s::world, "Smap"));
      (*Smap)["abij"] = (*Sm)["ai"] * (*Sm)["bj"];
      CTF::Bivar_Function<> fMultiply(&multiply<double>);
      Vpphh->contract(
        1.0, *Vpphh,"pqrs", *Smap,"pqrs", 0.0,"pqrs", fMultiply
      );
      delete Smap;
      delete Sm;
    }

    if ( isArgumentGiven("HHHHCoulombIntegrals") ) {
      allocatedTensorArgument<double>("HHHHCoulombIntegrals", Vhhhh);
    }
    allocatedTensorArgument<double>("PPHHCoulombIntegrals", Vpphh);

  }
  else {
    throw "Output either CoulombIntegrals or \
           PPHH && HHHH && OrbitalCoefficients && HoleEigenEnergies";
  }

}
