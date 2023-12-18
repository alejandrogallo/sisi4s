#include <algorithms/ParenthesisTriples.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/BinaryTensorFormat.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <fstream>
#include <extern/Lapack.hpp>
#include <math/ComplexTensor.hpp>
#include <numeric>
#include <util/Timer.hpp>
#include <algorithm>

using namespace sisi4s;

template <typename F>
IrmlerTensor<F> readBinaryTensorSerial(std::string filename,
                                       int64_t offset = 0,
                                       int64_t elements = 0,
                                       bool talk = true) {
  std::vector<F> data;
  std::vector<int64_t> lens;
  std::fstream f;
  int64_t length(1);
  BinaryTensorHeader header;
  if (talk)
    LOG(1, "SerialReader") << "reading tensor from  " << filename << std::endl;

  // read header to get order
  f.open(filename, std::ios::in | std::ios::binary);
  if (!f) throw EXCEPTION("file not found");
  f.read((char *)&header, sizeof(header));

  if (talk) LOG(1, "SerialReader") << "order " << header.order << std::endl;

  // for every dimension, read in the dimension header to find out the lens
  for (int i = 0; i < header.order; i++) {
    BinaryTensorDimensionHeader dimensionHeader;
    f.read((char *)&dimensionHeader, sizeof(dimensionHeader));
    if (talk)
      LOG(1, "SerialReader") << "lens " << dimensionHeader.length << std::endl;
    lens.push_back(dimensionHeader.length);
    length *= dimensionHeader.length;
  }

  if (offset > 0) {
    offset *= sizeof(F);
    offset += f.tellg(); // we have to add the offset of the header
    f.seekg(offset);
  }

  if (elements > 0) length = elements;

  // resize the data
  data.resize(length);
  f.read((char *)data.data(), length * sizeof(F));
  f.close();

  return IrmlerTensor<F>(lens, header.order, data);
}

// FUNCTIONS
double ParenthesisTriples::getEnergy(const double epsijk,
                                     const double *Tabc_,
                                     const double *Zabc_) {
  double energy(0.);
#pragma omp parallel for reduction(+ : energy)
  for (int64_t c = 0; c < Nv; c++)
    for (int64_t b = 0; b < Nv; b++) {
      double ec(epsa[c]);
      double eb(epsa[b]);
      for (int64_t a = 0; a < Nv; a++) {
        double ea(epsa[a]);
        double denominator(epsijk - ea - eb - ec);
        energy += Tabc_[a + Nv * b + Nv * Nv * c]
                * Zabc_[a + Nv * b + Nv * Nv * c] / denominator;
      }
    }
  return energy;
}

double getEnergyZero(const int Nv,
                     const double epsijk,
                     const double *epsa,
                     const double *Tabc_,
                     const double *Zabc_) {
  const int64_t blockSize = 16;
  double energy(0.);
#pragma omp parallel for reduction(+ : energy)
  for (int64_t cc = 0; cc < Nv; cc += blockSize) {
    int64_t cend = cc + blockSize < Nv ? blockSize : Nv - cc;
    for (int64_t bb(cc); bb < Nv; bb += blockSize) {
      int64_t bend = bb + blockSize < Nv ? blockSize : Nv - bb;
      for (int64_t aa(bb); aa < Nv; aa += blockSize) {
        int64_t aend = aa + blockSize < Nv ? blockSize : Nv - aa;
        for (int64_t c(cc); c < cc + cend; c++) {
          double ec(epsa[c]);
          int64_t bstart = bb > c ? bb : c;
          for (int64_t b(bstart); b < bb + bend; b++) {
            double eb(epsa[b]);
            double facbc(b == c ? 0.5 : 1.0);
            int64_t astart = aa > b ? aa : b;
            for (int64_t a(astart); a < aa + aend; a++) {
              double ea(epsa[a]);
              double facab(a == b ? 0.5 : 1.0);
              double denominator(epsijk - ea - eb - ec);
              double U(Zabc_[a + Nv * b + Nv * Nv * c]);
              double V(Zabc_[a + Nv * c + Nv * Nv * b]);
              double W(Zabc_[b + Nv * a + Nv * Nv * c]);
              double X(Zabc_[b + Nv * c + Nv * Nv * a]);
              double Y(Zabc_[c + Nv * a + Nv * Nv * b]);
              double Z(Zabc_[c + Nv * b + Nv * Nv * a]);

              double A(Tabc_[a + Nv * b + Nv * Nv * c]);
              double B(Tabc_[a + Nv * c + Nv * Nv * b]);
              double C(Tabc_[b + Nv * a + Nv * Nv * c]);
              double D(Tabc_[b + Nv * c + Nv * Nv * a]);
              double E(Tabc_[c + Nv * a + Nv * Nv * b]);
              double F(Tabc_[c + Nv * b + Nv * Nv * a]);
              double value(3.0 * (A * U + B * V + C * W + D * X + E * Y + F * Z)
                           + ((U + X + Y) - 2.0 * (V + W + Z)) * (A + D + E)
                           + ((V + W + Z) - 2.0 * (U + X + Y)) * (B + C + F));
              energy += 2.0 * value / denominator * facbc * facab;
            }
          }
        }
      }
    }
  }
  return energy;
}

double getEnergyOne(const int Nv,
                    const double epsijk,
                    const double *epsa,
                    const double *Tabc_,
                    const double *Zabc_) {
  const int64_t blockSize = 16;
  double energy(0.);
#pragma omp parallel for reduction(+ : energy)
  for (int64_t cc = 0; cc < Nv; cc += blockSize) {
    int64_t cend = cc + blockSize < Nv ? blockSize : Nv - cc;
    for (int64_t bb(cc); bb < Nv; bb += blockSize) {
      int64_t bend = bb + blockSize < Nv ? blockSize : Nv - bb;
      for (int64_t aa(bb); aa < Nv; aa += blockSize) {
        int64_t aend = aa + blockSize < Nv ? blockSize : Nv - aa;
        for (int64_t c(cc); c < cc + cend; c++) {
          double ec(epsa[c]);
          int64_t bstart = bb > c ? bb : c;
          for (int64_t b(bstart); b < bb + bend; b++) {
            double facbc(b == c ? 0.5 : 1.0);
            double eb(epsa[b]);
            int64_t astart = aa > b ? aa : b;
            for (int64_t a(astart); a < aa + aend; a++) {
              double ea(epsa[a]);
              double facab(a == b ? 0.5 : 1.0);
              double denominator(epsijk - ea - eb - ec);
              double U(Zabc_[a + Nv * b + Nv * Nv * c]);
              double V(Zabc_[b + Nv * c + Nv * Nv * a]);
              double W(Zabc_[c + Nv * a + Nv * Nv * b]);
              double A(Tabc_[a + Nv * b + Nv * Nv * c]);
              double B(Tabc_[b + Nv * c + Nv * Nv * a]);
              double C(Tabc_[c + Nv * a + Nv * Nv * b]);
              double value(3.0 * (A * U + B * V + C * W)
                           - (A + B + C) * (U + V + W));
              energy += 2.0 * value / denominator * facbc * facab;
            }
          }
        }
      }
    }
  }
  return energy;
}

inline void permuteAddOne(const int Nv, const double *input, double *output) {
#pragma omp parallel for
  for (int64_t k = 0; k < Nv; k++)
    for (int64_t j(0); j < Nv; j++)
      for (int64_t i(0); i < Nv; i++) {
        output[i + j * Nv + k * Nv * Nv] += input[i + j * Nv * Nv + k * Nv];
      }
}

inline void permuteMoveOne(const int Nv, const double *input, double *output) {
#pragma omp parallel for
  for (int64_t k = 0; k < Nv; k++)
    for (int64_t j(0); j < Nv; j++)
      for (int64_t i(0); i < Nv; i++) {
        output[i + j * Nv + k * Nv * Nv] = input[i + j * Nv * Nv + k * Nv];
      }
}

inline void permuteMoveTwo(const int Nv, const double *input, double *output) {
  const int64_t blockSize = 50;
#pragma omp parallel for
  for (int64_t k = (0); k < Nv; k++) {
    for (int64_t j(0); j < Nv; j += blockSize) {
      int64_t incj = j + blockSize < Nv ? blockSize : Nv - j;
      for (int64_t i(0); i < Nv; i += blockSize) {
        int64_t inci = i + blockSize < Nv ? blockSize : Nv - i;
        for (int64_t jj(j); jj < j + incj; jj++)
          for (int64_t ii(i); ii < i + inci; ii++) {
            output[ii + jj * Nv + k * Nv * Nv] =
                input[jj + ii * Nv + k * Nv * Nv];
          }
      }
    }
  }
}

inline void
fullPermutationZero(const int Nv, const double *input, double *output) {
  const int64_t blockSize = 16;
#pragma omp parallel for
  for (int64_t k = (0); k < Nv; k += blockSize) {
    int64_t inck = k + blockSize < Nv ? blockSize : Nv - k;
    for (int64_t j(0); j < Nv; j += blockSize) {
      int64_t incj = j + blockSize < Nv ? blockSize : Nv - j;
      for (int64_t i(0); i < Nv; i += blockSize) {
        int64_t inci = i + blockSize < Nv ? blockSize : Nv - i;
        for (int64_t kk(k); kk < k + inck; kk++)
          for (int64_t jj(j); jj < j + incj; jj++)
            for (int64_t ii(i); ii < i + inci; ii++) {
              output[ii + jj * Nv + kk * Nv * Nv] =
                  8.0 * input[ii + jj * Nv + kk * Nv * Nv]
                  - 4.0 * input[ii + kk * Nv + jj * Nv * Nv]
                  - 4.0 * input[jj + ii * Nv + kk * Nv * Nv]
                  + 2.0 * input[jj + kk * Nv + ii * Nv * Nv]
                  + 2.0 * input[kk + ii * Nv + jj * Nv * Nv]
                  - 4.0 * input[kk + jj * Nv + ii * Nv * Nv];
            }
      }
    }
  }
}

inline void
fullPermutationOne(const int Nv, const double *input, double *output) {
  const int64_t blockSize = 16;
#pragma omp parallel for
  for (int64_t k = (0); k < Nv; k += blockSize) {
    int64_t inck = k + blockSize < Nv ? blockSize : Nv - k;
    for (int64_t j(0); j < Nv; j += blockSize) {
      int64_t incj = j + blockSize < Nv ? blockSize : Nv - j;
      for (int64_t i(0); i < Nv; i += blockSize) {
        int64_t inci = i + blockSize < Nv ? blockSize : Nv - i;
        for (int64_t kk(k); kk < k + inck; kk++)
          for (int64_t jj(j); jj < j + incj; jj++)
            for (int64_t ii(i); ii < i + inci; ii++) {
              output[ii + jj * Nv + kk * Nv * Nv] =
                  2.0 * input[ii + jj * Nv + kk * Nv * Nv]
                  - 1.0 * input[ii + kk * Nv + jj * Nv * Nv]
                  - 1.0 * input[kk + ii * Nv + jj * Nv * Nv];
            }
      }
    }
  }
}

inline void
fullPermutationTwo(const int Nv, const double *input, double *output) {
  const int64_t blockSize = 16;
#pragma omp parallel for
  for (int64_t k = (0); k < Nv; k += blockSize) {
    int64_t inck = k + blockSize < Nv ? blockSize : Nv - k;
    for (int64_t j(0); j < Nv; j += blockSize) {
      int64_t incj = j + blockSize < Nv ? blockSize : Nv - j;
      for (int64_t i(0); i < Nv; i += blockSize) {
        int64_t inci = i + blockSize < Nv ? blockSize : Nv - i;
        for (int64_t kk(k); kk < k + inck; kk++)
          for (int64_t jj(j); jj < j + incj; jj++)
            for (int64_t ii(i); ii < i + inci; ii++) {
              output[ii + jj * Nv + kk * Nv * Nv] =
                  2.0 * input[ii + jj * Nv + kk * Nv * Nv]
                  - 1.0 * input[jj + ii * Nv + kk * Nv * Nv]
                  - 1.0 * input[jj + kk * Nv + ii * Nv * Nv];
            }
      }
    }
  }
}

void ParenthesisTriples::doublesContribution(const std::array<int64_t, 3> &ijk,
                                             IJKPointer integralContainer,
                                             double *scratchV,
                                             double *scratchO,
                                             double *output) {

  double mOne(-1.0);
  double pOne(1.0);
  double Zero(0.0);
  int64_t i(ijk[0]), j(ijk[1]), k(ijk[2]);
  // blas needs integer(maybe), these integer can not overflow for us:
  int blasNvNv(Nv * Nv), blasNv(Nv), blasNo(No);

  if (holeDiagram) {
// t[bali] * V[lcjk] = F[bac] -> T[abl](i) * V[lcjk] = F[abc]           H0
#pragma omp parallel for
    for (int64_t l = (0); l < No; l++) {
      for (int64_t b(0); b < Nv; b++) {
        for (int64_t a(0); a < Nv; a++) {
          scratchV[a + Nv * b + Nv * Nv * l] =
              Tabij[b + Nv * a + Nv * Nv * l + i * No * Nv * Nv];
        }
      }
    }
    dgemm_("N",
           "N",
           &blasNvNv,
           &blasNv,
           &blasNo,
           &mOne,
           scratchV,
           &blasNvNv,
           &Vhphh[j * No * Nv + k * No * No * Nv],
           &blasNo,
           &Zero,
           output,
           &blasNvNv);
    // t[ablj] * V[lcik] = F[abc]                                           H1
    dgemm_("N",
           "N",
           &blasNvNv,
           &blasNv,
           &blasNo,
           &mOne,
           &Tabij[j * No * Nv * Nv],
           &blasNvNv,
           &Vhphh[i * No * Nv + k * No * No * Nv],
           &blasNo,
           &pOne,
           output,
           &blasNvNv);
    // t[cali] * V[lbkj] = F[cab] -> T[acl](i)* V[lbkj] = F[acb] (PERM1)    H3
    dgemm_("N",
           "N",
           &blasNvNv,
           &blasNv,
           &blasNo,
           &mOne,
           scratchV,
           &blasNvNv,
           &Vhphh[k * No * Nv + j * No * No * Nv],
           &blasNo,
           &Zero,
           scratchO,
           &blasNvNv);
    // t[aclk] * V[lbij] = F[acb]                                (PERM1)   H5
    dgemm_("N",
           "N",
           &blasNvNv,
           &blasNv,
           &blasNo,
           &mOne,
           &Tabij[k * No * Nv * Nv],
           &blasNvNv,
           &Vhphh[i * No * Nv + j * No * No * Nv],
           &blasNo,
           &pOne,
           scratchO,
           &blasNvNv);
// t[cblj] * V[laki] = F[cba] -> V[laki] * T[bcl](j) = F[abc]           H4
#pragma omp parallel for
    for (int64_t l = (0); l < No; l++) {
      for (int64_t b(0); b < Nv; b++) {
        for (int64_t a(0); a < Nv; a++) {
          scratchV[a + Nv * b + Nv * Nv * l] =
              Tabij[b + Nv * a + Nv * Nv * l + j * No * Nv * Nv];
        }
      }
    }
    dgemm_("T",
           "T",
           &blasNv,
           &blasNvNv,
           &blasNo,
           &mOne,
           &Vhphh[k * No * Nv + i * No * No * Nv],
           &blasNo,
           scratchV,
           &blasNvNv,
           &pOne,
           output,
           &blasNv);
    // t[bclk] * V[laji] = F[bca] -> V[laji] * t[bclk] = F[abc]            H2
    dgemm_("T",
           "T",
           &blasNv,
           &blasNvNv,
           &blasNo,
           &mOne,
           &Vhphh[j * No * Nv + i * No * No * Nv],
           &blasNo,
           &Tabij[k * No * Nv * Nv],
           &blasNvNv,
           &pOne,
           output,
           &blasNv);
  }
  // PARTICLE CONTRACTION
  if (particleDiagram) {
    // when we dont calculate the holeDiagram we have to zero the
    // output/scratch0
    double zOrO = holeDiagram ? 1.0 : 0.0;
    // t[aeij] * V[bcek] = F[abc]                                          P0
    dgemm_("N",
           "T",
           &blasNv,
           &blasNvNv,
           &blasNv,
           &pOne,
           &Tabij[i * Nv * Nv + j * No * Nv * Nv],
           &blasNv,
           integralContainer.k,
           &blasNvNv,
           &zOrO,
           output,
           &blasNv);
    // t[ceki] * V[abej] = F[cab] -> V[abej] * t[ceki] = F[abc]            P5
    dgemm_("N",
           "T",
           &blasNvNv,
           &blasNv,
           &blasNv,
           &pOne,
           integralContainer.j,
           &blasNvNv,
           &Tabij[k * Nv * Nv + i * No * Nv * Nv],
           &blasNv,
           &pOne,
           output,
           &blasNvNv);
    // WE NEED A PERMUTED V[abe] (i)
    permuteMoveTwo(Nv, integralContainer.i, scratchV);
    // t[cekj] * V[baei] = F[cba] -> v[abe](i) * t[cekj] = F[abc]          P2
    dgemm_("N",
           "T",
           &blasNvNv,
           &blasNv,
           &blasNv,
           &pOne,
           scratchV,
           &blasNvNv,
           &Tabij[k * Nv * Nv + j * No * Nv * Nv],
           &blasNv,
           &pOne,
           output,
           &blasNvNv);
    // t[beji] * V[acek] = F[bac] -> V[acek] * t[beji] = F[acb] (PERM1)    P1
    dgemm_("N",
           "T",
           &blasNvNv,
           &blasNv,
           &blasNv,
           &pOne,
           integralContainer.k,
           &blasNvNv,
           &Tabij[j * Nv * Nv + i * No * Nv * Nv],
           &blasNv,
           &zOrO,
           scratchO,
           &blasNvNv);

    // t[aeik] * V[cbej] = F[acb]                               (PERM1)    P3
    dgemm_("N",
           "T",
           &blasNv,
           &blasNvNv,
           &blasNv,
           &pOne,
           &Tabij[i * Nv * Nv + k * No * Nv * Nv],
           &blasNv,
           integralContainer.j,
           &blasNvNv,
           &pOne,
           scratchO,
           &blasNv);
    // t[bejk] * V[caei] = F[bca] -> v[ace](i) * t[bejk] = F[acb] (PERM1)  P4
    dgemm_("N",
           "T",
           &blasNvNv,
           &blasNv,
           &blasNv,
           &pOne,
           scratchV,
           &blasNvNv,
           &Tabij[j * Nv * Nv + k * No * Nv * Nv],
           &blasNv,
           &pOne,
           scratchO,
           &blasNvNv);
  }
  permuteAddOne(Nv, scratchO, output);
}

void ParenthesisTriples::singlesContribution(const std::array<int64_t, 3> &ijk,
                                             double *scratchO,
                                             double *output) {
  int inc(1);
  int64_t i(ijk[0]), j(ijk[1]), k(ijk[2]);
  // dger expects integer arguments
  int blasNv(Nv), blasNvNv(Nv * Nv);
  double one(1.0);
  // ijk: T[ai]*V[bcjk] = F[abc]                                          S0
  // kij:  T[ak]*V[bcij] = F[abc] -> V[bcij]*T[ak]= F[bca]                ??
  // jik: T[aj]*V[bcik] = F[abc]
  int64_t offjk = Nv * Nv * j + Nv * Nv * No * k;
  int64_t offij = Nv * Nv * i + Nv * Nv * No * j;
  int64_t offik = Nv * Nv * i + Nv * Nv * No * k;

#pragma omp parallel for
  for (int64_t c = 0; c < Nv; c++)
    for (int64_t b = 0; b < Nv; b++)
      for (int64_t a = 0; a < Nv; a++) {
        output[a + b * Nv + c * Nv * Nv] +=
            Tai[a + i * Nv] * Vabij[offjk + b + c * Nv];
        output[a + b * Nv + c * Nv * Nv] +=
            Vabij[offij + a + b * Nv] * Tai[c + k * Nv];
        output[a + b * Nv + c * Nv * Nv] +=
            Vabij[offik + a + c * Nv] * Tai[b + j * Nv];
      }
}

void ParenthesisTriples::getVpppijkFromVertex(const std::array<int64_t, 3> &ijk,
                                              double *scratchO,
                                              double *Vpppijk) {
  double pOne(1.0), Zero(0.0);
  int64_t i(ijk[0]), j(ijk[1]), k(ijk[2]);
  // Target error [Nv*Nv, Nv], realGai is taken with the correct offset
  // ??? what does that comment mean??
  // Blas expects integer:
  int blasNvNv(Nv * Nv), blasNv(Nv), blasNG(NG);
  dgemm_("T",
         "N",
         &blasNvNv,
         &blasNv,
         &blasNG,
         &pOne,
         realGab,
         &blasNG,
         &realGai[i * Nv * NG],
         &blasNG,
         &Zero,
         scratchO,
         &blasNvNv);
  dgemm_("T",
         "N",
         &blasNvNv,
         &blasNv,
         &blasNG,
         &pOne,
         imagGab,
         &blasNG,
         &imagGai[i * Nv * NG],
         &blasNG,
         &pOne,
         scratchO,
         &blasNvNv);
  // We have to write the scratch on the container Vpppijk, however last two
  // indices flipped
  permuteMoveOne(Nv, scratchO, Vpppijk);

  dgemm_("T",
         "N",
         &blasNvNv,
         &blasNv,
         &blasNG,
         &pOne,
         realGab,
         &blasNG,
         &realGai[j * Nv * NG],
         &blasNG,
         &Zero,
         scratchO,
         &blasNvNv);
  dgemm_("T",
         "N",
         &blasNvNv,
         &blasNv,
         &blasNG,
         &pOne,
         imagGab,
         &blasNG,
         &imagGai[j * Nv * NG],
         &blasNG,
         &pOne,
         scratchO,
         &blasNvNv);
  // We have to write the scratch on the container Vpppijk, however last two
  // indices flipped
  permuteMoveOne(Nv, scratchO, &Vpppijk[NvCube]);

  dgemm_("T",
         "N",
         &blasNvNv,
         &blasNv,
         &blasNG,
         &pOne,
         realGab,
         &blasNG,
         &realGai[k * Nv * NG],
         &blasNG,
         &Zero,
         scratchO,
         &blasNvNv);
  dgemm_("T",
         "N",
         &blasNvNv,
         &blasNv,
         &blasNG,
         &pOne,
         imagGab,
         &blasNG,
         &imagGai[k * Nv * NG],
         &blasNG,
         &pOne,
         scratchO,
         &blasNvNv);
  // We have to write the scratch on the container Vpppijk, however last two
  // indices flipped
  permuteMoveOne(Nv, scratchO, &Vpppijk[2 * NvCube]);
}

IJKPointer
ParenthesisTriples::getVpppijkOnTheFly(const std::array<int64_t, 3> &ijk,
                                       VectorTensor<double> &vtensor) {
  // what do we need for ijk and its order
  std::list<int64_t> intersection;
  std::list<int64_t> difference;
  std::map<int64_t, IrmlerTensor<double> *> tensorMap;

  std::set_intersection(ijk.begin(),
                        ijk.end(),
                        vtensor.refs.begin(),
                        vtensor.refs.end(),
                        std::back_inserter(intersection));

  std::set_difference(vtensor.refs.begin(),
                      vtensor.refs.end(),
                      intersection.begin(),
                      intersection.end(),
                      std::back_inserter(difference));

  difference.unique();
  // only delete the tensor if it is not needed anyways
  // {006}->{011} have a difference of {06}
  for (const auto &ref : difference) {
    auto *p = std::find(ijk.begin(), ijk.end(), ref);
    if (std::distance(ijk.begin(), p) == 3) {
      delete vtensor.withReference(ref);
    }
  }

  for (const auto &ref : intersection) {
    tensorMap[ref] = vtensor.withReference(ref);
  }

  int64_t c = 0;
  for (const auto &i : ijk) {
    vtensor.refs[c] = i;
    if (tensorMap[i] == nullptr) {
      const int64_t offset = i * NvCube;
      vtensor.tensors[c] =
          new IrmlerTensor<double>(readBinaryTensorSerial<double>(
              in.get<std::string>("PPPHCoulombIntegralsFile",
                                  "PPPHCoulombIntegrals.bin"),
              offset,
              NvCube,
              false));
      tensorMap[i] = vtensor.withReference(i);
    } else {
      vtensor.tensors[c] = tensorMap[i];
    }
    c++;
  }

  return {vtensor.tensors[0]->data.data(),
          vtensor.tensors[1]->data.data(),
          vtensor.tensors[2]->data.data()};
}

DEFSPEC(
    ParenthesisTriples,
    SPEC_IN({"HoleEigenEnergiesFile",
             SPEC_VALUE_DEF("TODO: DOC", std::string, "HoleEigenEnergies.bin")},
            {"CoulombVertex",
             SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
    SPEC_OUT());

IMPLEMENT_ALGORITHM(ParenthesisTriples) {

  auto epsiT(readBinaryTensorSerial<double>(
      in.get<std::string>("HoleEigenEnergiesFile", "HoleEigenEnergies.bin")));
  epsi = epsiT.data.data();

  auto epsaT(readBinaryTensorSerial<double>(
      in.get<std::string>("ParticleEigenEnergiesFile",
                          "ParticleEigenEnergies.bin")));
  epsa = epsaT.data.data();

  No = epsiT.lens[0];
  Nv = epsaT.lens[0];

  NvCube = (int64_t)Nv * Nv * Nv;

  auto TaiT(readBinaryTensorSerial<double>(
      in.get<std::string>("CcsdSinglesAmplitudesFile",
                          "CcsdSinglesAmplitudes.bin")));
  Tai = TaiT.data.data();

  auto TabijT(readBinaryTensorSerial<double>(
      in.get<std::string>("CcsdDoublesAmplitudesFile",
                          "CcsdDoublesAmplitudes.bin")));
  Tabij = TabijT.data.data();

  auto VabijT(readBinaryTensorSerial<double>(
      in.get<std::string>("PPHHCoulombIntegralsFile",
                          "PPHHCoulombIntegrals.bin")));
  Vabij = VabijT.data.data();

  auto VijkaT(readBinaryTensorSerial<double>(
      in.get<std::string>("HPHHCoulombIntegralsFile",
                          "HPHHCoulombIntegrals.bin")));
  Vhphh = VijkaT.data.data();

  PTR(IrmlerTensor<double>) VppphT;

  VectorTensor<double> vtensor;

  bool fullPPPH(false);
  bool PPPHOnTheFly(false);
  if (isArgumentGiven("readPPPH")) {
    LOG(0, "Read and store PPPH CoulombIntegrals") << std::endl;
    VppphT = NEW(IrmlerTensor<double>,
                 readBinaryTensorSerial<double>(
                     in.get<std::string>("PPPHCoulombIntegralsFile",
                                         "PPPHCoulombIntegrals.bin")));
    Vppph = VppphT->data.data();
    fullPPPH = true;
  } else if (isArgumentGiven("PPPHOnTheFly")) {
    PPPHOnTheFly = true;
  } else {
    LOG(0, "Read and slice CoulombVertex") << std::endl;
    Tensor<sisi4s::complex> *GammaGqr(
        in.get<Tensor<sisi4s::complex> *>("CoulombVertex"));
    NG = GammaGqr->lens[0];
    int Np(GammaGqr->lens[1]);
    int iStart(0), iEnd(No);
    int aStart(Np - Nv), aEnd(Np);
    int GaiStart[] = {0, aStart, iStart};
    int GaiEnd[] = {NG, aEnd, iEnd};
    int GabStart[] = {0, aStart, aStart};
    int GabEnd[] = {NG, aEnd, aEnd};
    Tensor<sisi4s::complex> GammaGai(GammaGqr->slice(GaiStart, GaiEnd));
    Tensor<sisi4s::complex> GammaGab(GammaGqr->slice(GabStart, GabEnd));

    Tensor<double> realGammaGai(3,
                                GammaGai.lens,
                                GammaGai.sym,
                                *GammaGai.wrld,
                                "RealGammaGai");
    Tensor<double> imagGammaGai(3,
                                GammaGai.lens,
                                GammaGai.sym,
                                *GammaGai.wrld,
                                "ImagGammaGai");
    Tensor<double> realGammaGab(3,
                                GammaGab.lens,
                                GammaGab.sym,
                                *GammaGai.wrld,
                                "RealGammaGab");
    Tensor<double> imagGammaGab(3,
                                GammaGab.lens,
                                GammaGab.sym,
                                *GammaGab.wrld,
                                "ImagGammaGab");

    fromComplexTensor(GammaGai, realGammaGai, imagGammaGai);
    fromComplexTensor(GammaGab, realGammaGab, imagGammaGab);

    realGab = new double[NG * Nv * Nv];
    imagGab = new double[NG * Nv * Nv];
    realGai = new double[NG * Nv * No];
    imagGai = new double[NG * Nv * No];

    realGammaGai.read_all(realGai);
    imagGammaGai.read_all(imagGai);
    realGammaGab.read_all(realGab);
    imagGammaGab.read_all(imagGab);
  }

  if (isArgumentGiven("noParticleDiagram")) particleDiagram = false;
  if (isArgumentGiven("noHoleDiagram")) holeDiagram = false;

  // CONSTRUCT TUPLE LIST
  double energy(0.);
  // I=J=K is always zero. Number of non-zero tuples is given as below
  int64_t nTuples(No * (No + 1) * (No + 2) / 6 - No);
  // just for reducing panic: stdout every 10% of total work
  int64_t nMessages(9);
  int64_t iMessage(nTuples / nMessages);
  std::vector<int64_t> messages(nMessages);
  for (int64_t i(0); i < nMessages; i++) { messages[i] = iMessage * (i + 1); }
  std::reverse(messages.begin(), messages.end());

  // construct tuple list
  std::vector<std::array<int64_t, 3>> tuplesList(nTuples);
  int64_t u(0);
  for (int64_t i(0); i < No; i++)
    for (int64_t j(i); j < No; j++)
      for (int64_t k(j); k < No; k++) {
        if (i == j && i == k) continue;
        tuplesList[u][0] = i;
        tuplesList[u][1] = j;
        tuplesList[u][2] = k;
        u++;
      }

  // ALLOCATE THE THREE SCRATCHES WHICH WILL BE NEEDED
  Tabc = new double[NvCube]; // Doubles only (see piecuch)
  Zabc = new double[NvCube]; // Singles + Doubles (see piecuch)
  scratch = new double[NvCube];
  if (!fullPPPH && !PPPHOnTheFly) {
    Vpppijk = new double[3 * NvCube]; // stores the needed Vppph for given tuple
  }
  LOG(0, "START LOOP OVER") << u << " TUPLES " << std::endl;
  IJKPointer integralContainer;
  double VppphSeconds(0.0), doublesSeconds(0.0), permuteSeconds(0.0);
  double energySeconds(0.0), singlesSeconds(0.0);
  for (int64_t i = 0; i < nTuples; i++) {
    //  for (int64_t i = 0; i < 100; i++){
    // just provide some stdout
    if (i == messages.back()) {
      LOG(1, "Finished") << i << " out of total " << nTuples << " tuples"
                         << std::endl;
      messages.pop_back();
    }
    std::array<int64_t, 3> ijk;
    ijk[0] = tuplesList[i][0];
    ijk[1] = tuplesList[i][1];
    ijk[2] = tuplesList[i][2];
    int distinct(0);
    if (ijk[0] == ijk[1]) distinct++;
    if (ijk[1] == ijk[2]) distinct--;
    double epsijk(epsi[ijk[0]] + epsi[ijk[1]] + epsi[ijk[2]]);
    double tupleEnergy(0.);
    Time VppphTime;
    {
      Timer VppphTimer(&VppphTime);
      if (fullPPPH) {
        integralContainer.i = &Vppph[ijk[0] * NvCube];
        integralContainer.j = &Vppph[ijk[1] * NvCube];
        integralContainer.k = &Vppph[ijk[2] * NvCube];
        // #pragma omp parallel for
        //  for(int64_t v(0); v<NvCube; v++){
        //   Vpppijk = Vppph + ijk[0]*NvCube;
        //   Vpppijk + NvCube = Vppph + ijk[1]*NvCube;
        //   Vpppijk + NvCube + NvCube = Vppph + ijk[2]*NvCube;
        //   Vpppijk[v]          = Vppph[v + ijk[0]*NvCube];
        //   Vpppijk[v+NvCube]   = Vppph[v + ijk[1]*NvCube];
        //   Vpppijk[v+2*NvCube] = Vppph[v + ijk[2]*NvCube];
        // }
      } else if (PPPHOnTheFly) {
        integralContainer = getVpppijkOnTheFly(ijk, vtensor);
        // #pragma omp parallel for
        //  for (int64_t v(0); v<NvCube; v++){
        //   Vpppijk[v]          = integralContainer.i[v];
        //   Vpppijk[v+NvCube]   = integralContainer.j[v];
        //   Vpppijk[v+2*NvCube] = integralContainer.k[v];
        // }
      } else {
        getVpppijkFromVertex(ijk, scratch, Vpppijk);
        integralContainer.i = Vpppijk;
        integralContainer.j = &Vpppijk[NvCube];
        integralContainer.k = &Vpppijk[2 * NvCube];
      }
    }
    VppphSeconds += VppphTime.getFractionalSeconds();
    Time doublesTime;
    {
      Timer doublesTimer(&doublesTime);
      doublesContribution(ijk, integralContainer, Zabc, scratch, Tabc);
    }
    doublesSeconds += doublesTime.getFractionalSeconds();

    // we need the doubles-contributions two times. (s+d)*d

    Time singlesTime;
    {
      Timer singlesTimer(&singlesTime);
#pragma omp parallel for
      for (int64_t i = (0); i < NvCube; i++) Zabc[i] = Tabc[i];
      singlesContribution(ijk, scratch, Zabc);
    }
    singlesSeconds += singlesTime.getFractionalSeconds();
    Time energyTime;
    {
      Timer energyTimer(&energyTime);
      if (distinct == 0)
        tupleEnergy = getEnergyZero(Nv, epsijk, epsa, Tabc, Zabc);
      else tupleEnergy = getEnergyOne(Nv, epsijk, epsa, Tabc, Zabc);
    }
    energySeconds += energyTime.getFractionalSeconds();
    energy += tupleEnergy;
  }
  LOG(0, "LOOP FINISHED, energy") << energy << std::endl;
  LOG(0, "  Vppph Time") << VppphSeconds << std::endl;
  LOG(0, "doubles Time") << doublesSeconds << std::endl;
  LOG(0, "energy  Time") << energySeconds << std::endl;
  LOG(0, "singles Time") << singlesSeconds << std::endl;
  Scalar<double> ctfEnergy(*Sisi4s::world);
  ctfEnergy[""] = energy;
  // ctfEnergy.set_val(energy);  // ctfBug
  setRealArgument("Energy", ctfEnergy);

  // free up memory
  delete[] Tabc;
  delete[] Zabc;
  delete[] scratch;
  if (realGab) delete[] realGab;
  if (imagGab) delete[] imagGab;
  if (realGai) delete[] realGai;
  if (imagGai) delete[] imagGai;
  if (Vpppijk) delete[] Vpppijk;
}

void ParenthesisTriples::dryRun() {}
