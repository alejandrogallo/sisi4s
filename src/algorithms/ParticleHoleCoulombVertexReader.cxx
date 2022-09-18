#include <algorithms/ParticleHoleCoulombVertexReader.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <fstream>

using namespace sisi4s;

char const *ParticleHoleCoulombVertexReader::Header::MAGIC = "sisi4sFTOD";
char const *ParticleHoleCoulombVertexReader::Chunk::REALS_MAGIC = "FTODreal";
char const *ParticleHoleCoulombVertexReader::Chunk::IMAGS_MAGIC = "FTODimag";
char const *ParticleHoleCoulombVertexReader::Chunk::REALSIA_MAGIC = "FTIAreal";
char const *ParticleHoleCoulombVertexReader::Chunk::IMAGSIA_MAGIC = "FTIAimag";
char const *ParticleHoleCoulombVertexReader::Chunk::EPSILONS_MAGIC = "FTODepsi";

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombVertexReader);

ParticleHoleCoulombVertexReader::ParticleHoleCoulombVertexReader(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {

}

ParticleHoleCoulombVertexReader::~ParticleHoleCoulombVertexReader() {
}

struct Unrestricter {

  CTF::Tensor<sisi4s::complex>*
  doVertex(CTF::Tensor<sisi4s::complex> *GammaGqr) const {
    // The field variable NG remains the same
    int vertexLens[] = {static_cast<int>(GammaGqr->lens[0]),
                        static_cast<int>(2*GammaGqr->lens[1]),
                        static_cast<int>(2*GammaGqr->lens[2])};
    auto uGammaGqr(
      new Tensor<sisi4s::complex>(3, vertexLens, GammaGqr->sym, *Sisi4s::world, "uGammaGqr")
    );

    int *upUnrestrictedStates(new int[GammaGqr->lens[1]]);
    for (int q(0); q < GammaGqr->lens[1]; ++q) {
      upUnrestrictedStates[q] = 2*q;
    }
    int *upUp[] = { nullptr, upUnrestrictedStates, upUnrestrictedStates };
    // do uGammaGqr[G, upUn[q], upUn[r]] = GammaGqr[G,q,r] with upUn[q] = 2q.
    // NOTE: the behavior of all below permute calls is documented differently
    // in v1.4.1
    uGammaGqr->permute(1.0, *GammaGqr, upUp, 1.0);
    delete upUnrestrictedStates;

    int *downUnrestrictedStates(new int[GammaGqr->lens[1]]);
    for (int q(0); q < GammaGqr->lens[1]; ++q) {
      downUnrestrictedStates[q] = 2*q+1;
    }
    int *downDown[] = { nullptr, downUnrestrictedStates, downUnrestrictedStates };
    // do uGammaGqr[G, dnUn[q], dnUn[r]] = GammaGqr[G,q,r] with dnUn[q] = 2q+1.
    uGammaGqr->permute(1.0, *GammaGqr, downDown, 1.0);
    delete downUnrestrictedStates;

    // overwrite restricted vertex
    return uGammaGqr;
  }

  CTF::Tensor<double>*
  doEigenEnergies(CTF::Tensor<double> *eps) const {
    int lens[] = { static_cast<int>(2*eps->lens[0]) };
    auto uEps(
      new Tensor<double>(
        1,
        lens,
        eps->sym,
        *Sisi4s::world,
        ("u" + std::string{eps->get_name()}).c_str()
      )
    );

    int *upUnrestrictedStates(new int[eps->lens[0]]);
    for (int i(0); i < eps->lens[0]; ++i) {
      upUnrestrictedStates[i] = 2*i;
    }
    // do uEps[upUn[i]] = eps[i] with upUn[i] = 2i
    uEps->permute(1.0, *eps, &upUnrestrictedStates, 1.0);
    delete upUnrestrictedStates;

    int *downUnrestrictedStates(new int[eps->lens[0]]);
    for (int i(0); i < eps->lens[0]; ++i) {
      downUnrestrictedStates[i] = 2*i + 1;
    }
    // do uEps[dnUn[i]] = eps[i] with dnUn[i] = 2i+1
    uEps->permute(1.0, *eps, &downUnrestrictedStates, 1.0);
    delete downUnrestrictedStates;

    // overwrite restricted eigen energies
    // FIXME: eigen energies should be given after handleUnrestricted
    return uEps;
  }

};

void ParticleHoleCoulombVertexReader::run() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "Reader") <<
    "Reading Coulomb vertex from file " << fileName << std::endl;
  MPI_File file;
  int mpiError(
    MPI_File_open(
      Sisi4s::world->comm, fileName.c_str(), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &file
    )
  );
  if (mpiError) throw new EXCEPTION("Failed to open file");
  // Read header: obtain NG,No,Nv,Np
  Header header;
  MPI_Status status;
  MPI_File_read(file, &header, sizeof(header), MPI_BYTE, &status);
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new EXCEPTION("Invalid file format");
  int NG(header.NG);
  int No(header.No);
  int Nv(header.Nv);
  int Np(No + Nv);

  // Print NG, No, Nv, Np
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;

  // Allocate output tensors
  int vertexLens[] = { NG, Nv, No };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(No, *Sisi4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(Nv, *Sisi4s::world, "epsa"));
  Tensor<sisi4s::complex> *GammaGai(
    new Tensor<sisi4s::complex>(
      3, vertexLens, vertexSyms, *Sisi4s::world, "GammaGai"
    )
  );

  // Enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<sisi4s::complex>("ParticleHoleCoulombVertex", GammaGai);

  // Real and imaginary parts are read in seperately
  Tensor<> realGammaGai(
    3, vertexLens, vertexSyms, *Sisi4s::world, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
    3, vertexLens, vertexSyms, *Sisi4s::world, "ImagGammaGai"
  );

  int64_t offset(sizeof(header));
  MPI_Offset fileSize;
  MPI_File_get_size(file, &fileSize);
  Chunk chunk;
  while (offset < fileSize) {
    MPI_File_read_at(file, offset, &chunk, sizeof(chunk), MPI_BYTE, &status);
    if (strncmp(chunk.magic, Chunk::REALSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      // NOTE: dirty fix for bug in particle hole FTODDUMP
      chunk.size = sizeof(chunk) + sizeof(double) * Nv*No*NG;
      LOG(1, "Reader") << "reading " << realGammaGai.get_name() << std::endl;
      realGammaGai.read_dense_from_file(file, offset+sizeof(chunk));
    } else
    if (strncmp(chunk.magic, Chunk::IMAGSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      // NOTE: dirty fix for bug in particle hole FTODDUMP
      chunk.size = sizeof(chunk) + sizeof(double) * Nv*No*NG;
      LOG(1, "Reader") << "reading " << imagGammaGai.get_name() << std::endl;
      imagGammaGai.read_dense_from_file(file, offset+sizeof(chunk));
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(1, "Reader") << "reading " << epsi->get_name() << ", "
        << epsa->get_name() << std::endl;
      epsi->read_dense_from_file(file, offset+sizeof(chunk));
      epsa->read_dense_from_file(file, offset+sizeof(chunk)+No*sizeof(double));
    }
    offset += chunk.size;
  }
  MPI_File_close(&file);

  // Combine to complex tensor
  toComplexTensor(realGammaGai, imagGammaGai, *GammaGai);

  // handle unrestricted
  int unrestricted(getIntegerArgument("unrestricted", 0));
  Unrestricter u;
  if (unrestricted) {
    // Enter the allocated data (and by that type the output data to tensors)
    LOG(1, "Reader") << "Unrestricting " << epsi->get_name() << std::endl;
    allocatedTensorArgument("HoleEigenEnergies", u.doEigenEnergies(epsi));
    LOG(1, "Reader") << "Unrestricting " << epsa->get_name() << std::endl;
    allocatedTensorArgument("ParticleEigenEnergies", u.doEigenEnergies(epsa));
    LOG(1, "Reader") << "Unrestricting " << GammaGai->get_name() << std::endl;
    allocatedTensorArgument<complex>("ParticleHoleCoulombVertex", u.doVertex(GammaGai));
  }

}

void ParticleHoleCoulombVertexReader::dryRun() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "Reader") <<
    "Reading particle hole Coulomb vertex from file " << fileName << std::endl;
  std::ifstream file(fileName.c_str(), std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new EXCEPTION("Failed to open file");
  // Read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new EXCEPTION("Invalid file format");
  file.close();

  int NG(header.NG);
  int No(header.No);
  int Nv(header.Nv);
  int Np(No + Nv);

  // Print NG, No, Nv, Np
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;

  // Allocate output tensors
  int vertexLens[] = { NG, Nv, No };
  int vertexSyms[] = { NS, NS, NS };
  DryTensor<> *epsi(new DryVector<>(No, SOURCE_LOCATION));
  DryTensor<> *epsa(new DryVector<>(Nv, SOURCE_LOCATION));
  DryTensor<complex> *GammaGai(
    new DryTensor<complex>(3, vertexLens, vertexSyms, SOURCE_LOCATION)
  );
  // Enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ParticleHoleCoulombVertex", GammaGai
  );

  // Real and imaginary parts are read in seperately
  DryTensor<> realGammaGai(3, vertexLens, vertexSyms, SOURCE_LOCATION);
  DryTensor<> imagGammaGai(3, vertexLens, vertexSyms, SOURCE_LOCATION);
}
