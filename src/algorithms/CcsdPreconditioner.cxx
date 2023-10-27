#include <algorithms/CcsdPreconditioner.hpp>
#include <math/RandomTensor.hpp>
#include <util/MpiCommunicator.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/TensorIo.hpp>
#include <util/Exception.hpp>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/SharedPointer.hpp>

using namespace sisi4s;

template <typename F>
void filterOutSpinFlipEntries(Tensor<F> &t) {
  int64_t nValues;
  int64_t *globalIndices;
  F *values;
  int order(t.order);

  // read the values of the tensor in the current processor
  t.read_local(&nValues, &globalIndices, &values);

  for (int i = 0; i < nValues; i++) {
    int g = globalIndices[i];
    // global index "carry"
    int gc = g;
    // this is a selector for the case where crossterms appear
    bool up(true), down(true);
    for (int o(0); o < order; o++) {
      int modulizer = t.lens[o];
      int ijk = gc % modulizer;
      up = up && (ijk % 2 == 0);
      down = down && ((ijk % 2 + 1) == 1);
      gc /= modulizer;
    }
    // if up and down are zero, it means that at some point the the indices
    // were part of a cross-term element, this should the be set to zero
    // because we're filtering out the terms of this kind.
    if (!(up || down)) { values[i] = F(0); }
  }

  // now we have to write the values back into the tensor
  t.write(nValues, globalIndices, values);

  // clean up the mess
  delete[] values;
  delete[] globalIndices;
}

/**
 * \brief Comparator that should filter out zero values of the diagonal
 * matrix.
 * Zero values are treated as infinite so that they get appended to the
 * end of the list.
 */
template <typename F>
class EomDiagonalValueComparator {
public:
  bool operator()(const std::pair<int, F> &a, const std::pair<int, F> &b) {
    F A(std::abs(a.second) < 1E-13 ? std::numeric_limits<F>::infinity()
                                   : a.second);
    F B(std::abs(b.second) < 1E-13 ? std::numeric_limits<F>::infinity()
                                   : b.second);
    double diff(computeDifference(A, B));
    // maintain magnitude finite!
    double magnitude(std::abs(a.second) + std::abs(b.second));
    if (diff > +1E-13 * magnitude) return true;
    if (diff < -1E-13 * magnitude) return false;
    return a.first < b.first;
  }

  double computeDifference(const F &a, const F &b) { return b - a; }
};

template <>
double EomDiagonalValueComparator<sisi4s::complex>::computeDifference(
    const sisi4s::complex &a,
    const sisi4s::complex &b) {
  double diff(b.imag() + b.real() - a.imag() - a.real());
  return diff;
}

template <typename F>
void CcsdPreconditioner<F>::calculateDiagonal() {
  diagonalH = NEW(V,
                  std::vector<PTR(Tensor<F>)>(
                      {NEW(Tensor<F>, *Tai), NEW(Tensor<F>, *Tabij)}),
                  std::vector<std::string>({"ai", "abij"}));
  // pointers to singles and doubles tensors of diagonal part
  auto Dai(diagonalH->get(0));
  auto Dabij(diagonalH->get(1));

  // TODO: Maybe insert the Tai part to the diagonal

  // calculate diagonal elements of H
  (*Dai)["bi"] = (-1.0) * (*Fij)["ii"];
  (*Dai)["bi"] += (+1.0) * (*Fab)["bb"];

  /*
    if (Viajb) {
      (*Dai)["bi"] += ( - 1.0 ) * (*Viajb)["ibib"];
    }
    if (Vijab) {
      (*Dai)["bi"] += ( + 1.0 ) * (*Tabij)["cbli"] * (*Vijab)["licb"];
      (*Dai)["bi"] += ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["micd"];
      (*Dai)["bi"] += ( - 0.5 ) * (*Tabij)["cblm"] * (*Vijab)["lmcb"];
    }
  */
  (*Dabij)["cdij"] = (-1.0) * (*Fij)["ii"];
  (*Dabij)["cdij"] += (-1.0) * (*Fij)["jj"];
  (*Dabij)["cdij"] += (+1.0) * (*Fab)["cc"];
  (*Dabij)["cdij"] += (+1.0) * (*Fab)["dd"];
  /*
    if (Vijkl) (*Dabij)["cdij"] += ( + 0.5 ) * (*Vijkl)["ijij"];

    if (Viajb) {
      (*Dabij)["ccij"] += ( + 1.0 ) * (*Viajb)["icic"];
      (*Dabij)["cdij"] += ( - 1.0 ) * (*Viajb)["icic"];
      (*Dabij)["ccii"] += ( - 1.0 ) * (*Viajb)["icic"];
      (*Dabij)["cdii"] += ( + 1.0 ) * (*Viajb)["icic"];
    }

    if (Vabcd) (*Dabij)["cdij"] += ( + 0.5 ) * (*Vabcd)["cdcd"];

    if (Vijab) {
      (*Dabij)["ccij"] += ( + 0.5 ) * (*Tabij)["ecij"] * (*Vijab)["ijec"];
      (*Dabij)["cdij"] += ( - 0.5 ) * (*Tabij)["ecij"] * (*Vijab)["ijec"];
      (*Dabij)["cdij"] += ( + 0.25) * (*Tabij)["efij"] * (*Vijab)["ijef"];
      (*Dabij)["cdij"] += ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["micd"];
      (*Dabij)["cdii"] += ( + 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["micd"];
      (*Dabij)["ccij"] += ( - 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["niec"];
      (*Dabij)["cdij"] += ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["niec"];
      (*Dabij)["ccii"] += ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["niec"];
      (*Dabij)["cdii"] += ( - 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["niec"];
      (*Dabij)["cdij"] += ( - 0.5 ) * (*Tabij)["efoi"] * (*Vijab)["oief"];
      (*Dabij)["cdii"] += ( + 0.5 ) * (*Tabij)["efoi"] * (*Vijab)["oief"];
      (*Dabij)["cdij"] += ( + 0.25) * (*Tabij)["cdmn"] * (*Vijab)["mncd"];
      (*Dabij)["ccij"] += ( + 0.5 ) * (*Tabij)["ecno"] * (*Vijab)["noec"];
      (*Dabij)["cdij"] += ( - 0.5 ) * (*Tabij)["ecno"] * (*Vijab)["noec"];
    }
  */
}

template <typename F>
std::vector<SDFockVector<F>>
CcsdPreconditioner<F>::getInitialBasis(const int eigenVectorsCount) {
  calculateDiagonal();
  LOG(0, "CcsdPreconditioner") << "Getting initial basis " << std::endl;
  if (preconditionerRandom) {
    LOG(0, "CcsdPreconditioner") << "Randomizing initial guess" << std::endl;
  }
  DefaultRandomEngine randomEngine;
  std::normal_distribution<double> normalDistribution(
      0.0,
      preconditionerRandomSigma);
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements(diagonalH->readLocal());
  // sort the local elements according to the eom comparator
  std::sort(localElements.begin(),
            localElements.end(),
            EomDiagonalValueComparator<F>());
  int localElementsSize(localElements.size());

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10 * eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  // get the local elements indices and values into their own vectors
  for (int i(0); i < std::min(localElementsSize, trialEigenVectorsCount); ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Sisi4s::world);
  std::vector<size_t> lowestElementIndices;
  std::vector<F> lowestElementValues;
  // get local lowest (indices or values) into a vector
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementValues.size());
  for (unsigned int i(0); i < lowestElementValues.size(); ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(lowestElements.begin(),
            lowestElements.end(),
            EomDiagonalValueComparator<F>());
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;

  int currentEigenVectorCount(0);
  unsigned int loopCount(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    V basisElement(*diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t, F>> elements;
    if (communicator.getRank() == 0) {
      if (loopCount >= lowestElements.size()) {
        throw EXCEPTION("No more elements to create initial basis");
      }
      // put a a pair in the elements (globalIndex, 1.0)
      elements.push_back(std::make_pair(lowestElements[loopCount].first, 1.0));
    }
    basisElement.write(elements);

    // random transformation
    if (preconditionerRandom) {
      auto Rai(*basisElement.get(0));
      auto Rabij(*basisElement.get(1));
      setRandomTensor(Rai, normalDistribution, randomEngine);
      setRandomTensor(Rabij, normalDistribution, randomEngine);
      (*basisElement.get(0))["ai"] += Rai["ai"];
      (*basisElement.get(1))["abij"] += Rabij["abij"];
    }

    // FILTER: unphysical components from the basisElement
    (*basisElement.get(1))["abii"] = 0.0;
    (*basisElement.get(1))["aaij"] = 0.0;
    (*basisElement.get(1))["aaii"] = 0.0;

    // FILTER: Antisymmetrize the new basis element
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["abji"];
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["baij"];

    // FILTER: Spin Flip
    if (!spinFlip) {
      LOG(0, "CcsdPreconditioner") << "Filtering out spin flip" << std::endl;
      for (auto &t : basisElement.component_tensors) {
        filterOutSpinFlipEntries(*t);
      }
    }

    // FILTER: Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basis[j].dot(basisElement);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if (std::abs(basisElementNorm) < 1e-10) {
      zeroVectorCount++;
      loopCount++;
      continue;
    }

    OUT() << "\tbasisSize=" << basis.size() << std::endl;

    basisElement =
        1.0 / std::sqrt(basisElement.dot(basisElement)) * basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    loopCount++;

    if (std::abs(basisElementNorm - double(1)) > 1e-10 * double(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);
  }

  // Now make sure that you are giving an antisymmetric basis
  // and also that it is grammschmited afterwards
  // LOG(1,"CcsdPreconditioner") << "Antisymmetrize basis" << std::endl;
  // for (unsigned int j(0); j < basis.size(); ++j) {
  //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["abji"];
  //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["baij"];
  //}

  // LOG(1,"CcsdPreconditioner") <<
  //"Performing Gramm Schmidt in the initial basis"
  //<< std::endl;
  // for (unsigned int b(0); b < basis.size(); ++b) {
  // V newVector(basis[b]);
  // for (unsigned int j(0); j < b; ++j) {
  // newVector -= basis[j] * basis[j].dot(basis[b]);
  //}
  //// normalize
  // basis[b] = 1 / std::sqrt(newVector.dot(newVector)) * newVector;
  //}

  OUT() << "\treturning vecs =" << basis.size() << std::endl;

  return basis;
}

template <typename F>
SDTFockVector<F>
CcsdPreconditioner<F>::getCorrection(const sisi4s::complex lambda,
                                     SDTFockVector<F> &residuum) {
  // Cast ccsdt into ccsd
  V w(residuum);
  // apply the old getCorrection
  // and cast again into a ccsdt vector
  SDTFockVector<F> result(getCorrection(lambda, w));
  // et voila
  return result;
}

template <typename F>
SFockVector<F>
CcsdPreconditioner<F>::getCorrection(const sisi4s::complex lambda,
                                     SFockVector<F> &residuum) {
  // Cast s into sd
  SDFockVector<F> w(residuum);
  // apply the old getCorrection
  // and cast again into a ccsdt vector
  SFockVector<F> result(getCorrection(lambda, w));
  // et voila
  return result;
}

template <typename F>
SDFockVector<F>
CcsdPreconditioner<F>::getCorrection(const sisi4s::complex lambda,
                                     SDFockVector<F> &residuum) {
  SDFockVector<F> w(*diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
  public:
    DiagonalCorrection(const double lambda_)
        : lambda(lambda_) {}
    F operator()(const F residuumElement, const F diagonalElement) {
      return std::abs(lambda - diagonalElement) < 1E-4
               ? 0.0
               : residuumElement / (lambda - diagonalElement);
    }

  protected:
    double lambda;
  } diagonalCorrection(std::real(lambda));

  SDFockVector<F> correction(*diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.get_components_count(); ++c) {
    const char *indices(correction.component_indices[c].c_str());
    (*correction.get(c))
        .contract(1.0,
                  *residuum.get(c),
                  indices,
                  *diagonalH->get(c),
                  indices,
                  0.0,
                  indices,
                  CTF::Bivar_Function<F>(diagonalCorrection));
  }
  // Filter out unphysical components from the correction
  (*correction.get(1))["abii"] = 0.0;
  (*correction.get(1))["aaij"] = 0.0;
  (*correction.get(1))["aaii"] = 0.0;

  // Antisymmetrize the correction
  (*correction.get(1))["abij"] -= (*correction.get(1))["abji"];
  (*correction.get(1))["abij"] -= (*correction.get(1))["baij"];
  (*correction.get(1))["abij"] = 0.25 * (*correction.get(1))["abij"];

  return correction;
}

template <typename F>
void EACcsdPreconditioner<F>::calculateDiagonal() {
  int No(this->Fij->lens[0]), Nv(this->Fab->lens[0]);
  std::vector<int> v{{Nv}}, vvo{{Nv, Nv, No}}, ns{{NS}}, nss{{NS, NS, NS}};

  auto &Fij = this->Fij;
  auto &Fab = this->Fab;

  this->diagonalH = NEW(
      SDFockVector<F>,
      std::vector<PTR(Tensor<F>)>(
          {NEW(Tensor<F>, 1, v.data(), ns.data(), *Sisi4s::world, "Da"),
           NEW(Tensor<F>, 3, vvo.data(), nss.data(), *Sisi4s::world, "Dabi")}),
      std::vector<std::string>({"a", "abi"}));

  auto Da(this->diagonalH->get(0));
  auto Dabi(this->diagonalH->get(1));

  // calculate diagonal elements of H
  (*Da)["a"] = (+1.0) * (*Fab)["aa"];

  (*Dabi)["cdi"] = (-1.0) * (*Fij)["ii"];
  (*Dabi)["cdi"] += (-1.0) * (*Fij)["jj"];
  (*Dabi)["cdi"] += (+1.0) * (*Fab)["cc"];
  (*Dabi)["cdi"] += (+1.0) * (*Fab)["dd"];
}

template <typename F>
std::vector<SDFockVector<F>>
EACcsdPreconditioner<F>::getInitialBasis(const int eigenVectorsCount) {
  calculateDiagonal();
  LOG(0, "EACcsdPreconditioner") << "Getting initial basis " << std::endl;
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements(this->diagonalH->readLocal());
  std::sort(localElements.begin(),
            localElements.end(),
            EomDiagonalValueComparator<F>());
  int localElementsSize(localElements.size());

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10 * eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  for (int i(0); i < std::min(localElementsSize, trialEigenVectorsCount); ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Sisi4s::world);
  std::vector<size_t> lowestElementIndices;
  std::vector<F> lowestElementValues;
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementValues.size());
  for (unsigned int i(0); i < lowestElementValues.size(); ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(lowestElements.begin(),
            lowestElements.end(),
            EomDiagonalValueComparator<F>());
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<SDFockVector<F>> basis;

  int currentEigenVectorCount(0);
  unsigned int b(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    SDFockVector<F> basisElement(*this->diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t, F>> elements;
    if (communicator.getRank() == 0) {
      if (b >= lowestElements.size()) {
        throw EXCEPTION("No more elements to create initial basis");
      }
      elements.push_back(std::make_pair(lowestElements[b].first, 1.0));
    }
    basisElement.write(elements);
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    (*basisElement.get(1))["aai"] = 0.0;

    // Antisymmetrize the new basis element
    (*basisElement.get(1))["abi"] -= (*basisElement.get(1))["bai"];

    LOG(1, "EACcsdPreconditioner")
        << "basis size " << basis.size() << std::endl;

    // Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basis[j].dot(basisElement);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if (std::abs(basisElementNorm) < 1e-10) {
      zeroVectorCount++;
      b++;
      continue;
    }

    basisElement =
        1.0 / std::sqrt(basisElement.dot(basisElement)) * basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    b++;

    if (std::abs(basisElementNorm - double(1)) > 1e-10 * double(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);
  }

  return basis;
}

template <typename F>
SDFockVector<F>
EACcsdPreconditioner<F>::getCorrection(const sisi4s::complex lambda,
                                       SDFockVector<F> &residuum) {
  SDFockVector<F> w(*this->diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
  public:
    DiagonalCorrection(const double lambda_)
        : lambda(lambda_) {}
    F operator()(const F residuumElement, const F diagonalElement) {
      return std::abs(lambda - diagonalElement) < 1E-4
               ? 0.0
               : residuumElement / (lambda - diagonalElement);
    }

  protected:
    double lambda;
  } diagonalCorrection(std::real(lambda));

  SDFockVector<F> correction(*this->diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.get_components_count(); ++c) {
    const char *indices(correction.component_indices[c].c_str());
    (*correction.get(c))
        .contract(1.0,
                  *residuum.get(c),
                  indices,
                  *this->diagonalH->get(c),
                  indices,
                  0.0,
                  indices,
                  CTF::Bivar_Function<F>(diagonalCorrection));
  }
  // Filter out unphysical components from the correction
  (*correction.get(1))["aai"] = 0.0;

  // Antisymmetrize the correction
  (*correction.get(1))["abi"] -= (*correction.get(1))["bai"];

  return correction;
}

template <typename F>
void IPCcsdPreconditioner<F>::calculateDiagonal() {
  int No(this->Fij->lens[0]), Nv(this->Fab->lens[0]);
  std::vector<int> o{{No}}, voo{{Nv, No, No}}, ns{{NS}}, nss{{NS, NS, NS}};

  auto &Fij = this->Fij;
  auto &Fab = this->Fab;

  this->diagonalH = NEW(
      SDFockVector<F>,
      std::vector<PTR(Tensor<F>)>(
          {NEW(Tensor<F>, 1, o.data(), ns.data(), *Sisi4s::world, "Di"),
           NEW(Tensor<F>, 3, voo.data(), nss.data(), *Sisi4s::world, "Daij")}),
      std::vector<std::string>({"i", "aij"}));

  auto Di(this->diagonalH->get(0));
  auto Daij(this->diagonalH->get(1));

  // calculate diagonal elements of H
  (*Di)["i"] = (-1.0) * (*Fij)["ii"];

  (*Daij)["cij"] = (-1.0) * (*Fij)["ii"];
  (*Daij)["cij"] += (-1.0) * (*Fij)["jj"];
  (*Daij)["cij"] += (+1.0) * (*Fab)["cc"];
  (*Daij)["cij"] += (+1.0) * (*Fab)["dd"];
}

template <typename F>
std::vector<SDFockVector<F>>
IPCcsdPreconditioner<F>::getInitialBasis(const int eigenVectorsCount) {
  calculateDiagonal();
  LOG(0, "IPCcsdPreconditioner") << "Getting initial basis " << std::endl;
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements(this->diagonalH->readLocal());
  std::sort(localElements.begin(),
            localElements.end(),
            EomDiagonalValueComparator<F>());
  int localElementsSize(localElements.size());

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10 * eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  for (int i(0); i < std::min(localElementsSize, trialEigenVectorsCount); ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Sisi4s::world);
  std::vector<size_t> lowestElementIndices;
  std::vector<F> lowestElementValues;
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementValues.size());
  for (unsigned int i(0); i < lowestElementValues.size(); ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(lowestElements.begin(),
            lowestElements.end(),
            EomDiagonalValueComparator<F>());
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<SDFockVector<F>> basis;

  int currentEigenVectorCount(0);
  unsigned int b(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    SDFockVector<F> basisElement(*this->diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t, F>> elements;
    if (communicator.getRank() == 0) {
      if (b >= lowestElements.size()) {
        throw EXCEPTION("No more elements to create initial basis");
      }
      elements.push_back(std::make_pair(lowestElements[b].first, 1.0));
    }
    basisElement.write(elements);
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    (*basisElement.get(1))["aii"] = 0.0;

    // Antisymmetrize the new basis element
    (*basisElement.get(1))["aij"] -= (*basisElement.get(1))["aji"];

    LOG(1, "IPCcsdPreconditioner")
        << "basis size " << basis.size() << std::endl;

    // Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basis[j].dot(basisElement);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if (std::abs(basisElementNorm) < 1e-10) {
      zeroVectorCount++;
      b++;
      continue;
    }

    basisElement =
        1.0 / std::sqrt(basisElement.dot(basisElement)) * basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    b++;

    if (std::abs(basisElementNorm - double(1)) > 1e-10 * double(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);
  }

  return basis;
}

template <typename F>
SDFockVector<F>
IPCcsdPreconditioner<F>::getCorrection(const sisi4s::complex lambda,
                                       SDFockVector<F> &residuum) {
  SDFockVector<F> w(*this->diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
  public:
    DiagonalCorrection(const double lambda_)
        : lambda(lambda_) {}
    F operator()(const F residuumElement, const F diagonalElement) {
      return std::abs(lambda - diagonalElement) < 1E-4
               ? 0.0
               : residuumElement / (lambda - diagonalElement);
    }

  protected:
    double lambda;
  } diagonalCorrection(std::real(lambda));

  SDFockVector<F> correction(*this->diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.get_components_count(); ++c) {
    const char *indices(correction.component_indices[c].c_str());
    (*correction.get(c))
        .contract(1.0,
                  *residuum.get(c),
                  indices,
                  *this->diagonalH->get(c),
                  indices,
                  0.0,
                  indices,
                  CTF::Bivar_Function<F>(diagonalCorrection));
  }
  // Filter out unphysical components from the correction
  (*correction.get(1))["aii"] = 0.0;

  // Antisymmetrize the correction
  (*correction.get(1))["aij"] -= (*correction.get(1))["aji"];

  return correction;
}

template <typename F>
void CISPreconditioner<F>::calculateDiagonal() {
  int No(this->Fij->lens[0]), Nv(this->Fab->lens[0]);
  std::vector<int> vo{{Nv, No}}, ns{{NS, NS}};

  auto &Fij = this->Fij;
  auto &Fab = this->Fab;

  this->diagonalH =
      NEW(typename CISPreconditioner<F>::V,
          std::vector<PTR(Tensor<F>)>(
              {NEW(Tensor<F>, 2, vo.data(), ns.data(), *Sisi4s::world, "Dai")}),
          std::vector<std::string>({"ai"}));

  auto Dai(this->diagonalH->get(0));
  // auto Dabij(this->diagonalH->get(1));

  // calculate diagonal elements of H
  (*Dai)["bi"] = (-1.0) * (*Fij)["ii"];
  (*Dai)["bi"] += (+1.0) * (*Fab)["bb"];
  //   (*Dabij)["cdij"] = (-1.0) * (*Fij)["ii"];
  //   (*Dabij)["cdij"] += (-1.0) * (*Fij)["jj"];
  //   (*Dabij)["cdij"] += (+1.0) * (*Fab)["cc"];
  //   (*Dabij)["cdij"] += (+1.0) * (*Fab)["dd"];
}

template <typename F>
typename CISPreconditioner<F>::V CISPreconditioner<F>::getCorrection(
    const sisi4s::complex lambda,
    typename CISPreconditioner<F>::V &residuum) {
  typename CISPreconditioner<F>::V w(*this->diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
  public:
    DiagonalCorrection(const double lambda_)
        : lambda(lambda_) {}
    F operator()(const F residuumElement, const F diagonalElement) {
      return std::abs(lambda - diagonalElement) < 1E-4
               ? 0.0
               : residuumElement / (lambda - diagonalElement);
    }

  protected:
    double lambda;
  } diagonalCorrection(std::real(lambda));

  typename CISPreconditioner<F>::V correction(*this->diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.get_components_count(); ++c) {
    const char *indices(correction.component_indices[c].c_str());
    (*correction.get(c))
        .contract(1.0,
                  *residuum.get(c),
                  indices,
                  *this->diagonalH->get(c),
                  indices,
                  0.0,
                  indices,
                  CTF::Bivar_Function<F>(diagonalCorrection));
  }
  return correction;
}

template <typename F>
std::vector<typename CISPreconditioner<F>::V>
CISPreconditioner<F>::getInitialBasis(const int eigenVectorsCount) {
  calculateDiagonal();
  LOG(0, "CISPreconditioner") << "Getting initial basis " << std::endl;
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements(this->diagonalH->readLocal());
  std::sort(localElements.begin(),
            localElements.end(),
            EomDiagonalValueComparator<F>());
  int localElementsSize(localElements.size());

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10 * eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  for (int i(0); i < std::min(localElementsSize, trialEigenVectorsCount); ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Sisi4s::world);
  std::vector<size_t> lowestElementIndices;
  std::vector<F> lowestElementValues;
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementValues.size());
  for (unsigned int i(0); i < lowestElementValues.size(); ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(lowestElements.begin(),
            lowestElements.end(),
            EomDiagonalValueComparator<F>());
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<typename CISPreconditioner<F>::V> basis;

  int currentEigenVectorCount(0);
  unsigned int b(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    typename CISPreconditioner<F>::V basisElement(*this->diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t, F>> elements;
    if (communicator.getRank() == 0) {
      if (b >= lowestElements.size()) {
        throw EXCEPTION("No more elements to create initial basis");
      }
      elements.push_back(std::make_pair(lowestElements[b].first, 1.0));
    }
    basisElement.write(elements);
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    // (*basisElement.get(1))["aii"] = 0.0;

    // Antisymmetrize the new basis element
    // (*basisElement.get(1))["aij"] -= (*basisElement.get(1))["aji"];

    LOG(1, "CISPreconditioner") << "basis size " << basis.size() << std::endl;

    // Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basis[j].dot(basisElement);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if (std::abs(basisElementNorm) < 1e-10) {
      zeroVectorCount++;
      b++;
      continue;
    }

    basisElement =
        1.0 / std::sqrt(basisElement.dot(basisElement)) * basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    b++;

    if (std::abs(basisElementNorm - double(1)) > 1e-10 * double(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);
  }

  return basis;
}

// instantiate
template class sisi4s::CcsdPreconditioner<double>;
template class sisi4s::CcsdPreconditioner<sisi4s::complex>;

template class sisi4s::IPCcsdPreconditioner<double>;
template class sisi4s::IPCcsdPreconditioner<sisi4s::complex>;

template class sisi4s::EACcsdPreconditioner<double>;
template class sisi4s::EACcsdPreconditioner<sisi4s::complex>;

template class sisi4s::CISPreconditioner<double>;
template class sisi4s::CISPreconditioner<sisi4s::complex>;
