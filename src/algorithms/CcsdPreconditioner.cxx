#include <algorithms/CcsdPreconditioner.hpp>
#include <math/RandomTensor.hpp>
#include <util/MpiCommunicator.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/TensorIo.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>


using namespace cc4s;

/**
 * \brief Comparator that should filter out zero values of the diagonal
 * matrix.
 * Zero values are treated as infinite so that they get appended to the
 * end of the list.
 */
template <typename F>
class EomDiagonalValueComparator {
public:
  bool operator ()(
    const std::pair<int, F> &a,
    const std::pair<int, F> &b
  ) {
    F A(
      std::abs(a.second) < 1E-13 ?
        std::numeric_limits<F>::infinity() : a.second
    );
    F B(
      std::abs(b.second) < 1E-13 ?
        std::numeric_limits<F>::infinity() : b.second
    );
    double diff(computeDifference(A, B));
    // maintain magnitude finite!
    double magnitude(std::abs(a.second)+std::abs(b.second));
    if (diff > +1E-13*magnitude) return true;
    if (diff < -1E-13*magnitude) return false;
    return a.first < b.first;
  }

  double computeDifference(const F &a, const F &b) { return b - a; }

};

template<>
double EomDiagonalValueComparator<complex>::computeDifference(
    const complex &a,
    const complex &b
  ) {
  double diff(b.imag() + b.real() - a.imag() - a.real());
  return diff;
}

template <typename F>
CcsdPreconditioner<F>::CcsdPreconditioner(
  CTF::Tensor<F> &Tai,
  CTF::Tensor<F> &Tabij,
  CTF::Tensor<F> &Fij,
  CTF::Tensor<F> &Fab,
  CTF::Tensor<F> &Vabcd,
  CTF::Tensor<F> &Viajb,
  CTF::Tensor<F> &Vijab,
  CTF::Tensor<F> &Vijkl
): diagonalH(
    std::vector<PTR(CTF::Tensor<F>)>(
      {NEW(CTF::Tensor<F>, Tai), NEW(CTF::Tensor<F>, Tabij)}
    ),
    std::vector<std::string>({"ai", "abij"})
  )
  {
  // pointers to singles and doubles tensors of diagonal part
  auto Dai( diagonalH.get(0) );
  auto Dabij( diagonalH.get(1) );

  // TODO: Maybe inster the Tai part to the diagonal

  // calculate diagonal elements of H
  (*Dai)["bi"] =  ( - 1.0 ) * Fij["ii"];
  (*Dai)["bi"] += ( + 1.0 ) * Fab["bb"];
/*  (*Dai)["bi"] += ( - 1.0 ) * Viajb["ibib"];
  (*Dai)["bi"] += ( + 1.0 ) * Tabij["cbli"] * Vijab["licb"];
  (*Dai)["bi"] += ( - 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dai)["bi"] += ( - 0.5 ) * Tabij["cblm"] * Vijab["lmcb"];
*/
  (*Dabij)["cdij"] =  ( - 1.0 ) * Fij["ii"];
  (*Dabij)["cdij"] += ( - 1.0 ) * Fij["jj"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Fab["cc"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Fab["dd"];
/*
  (*Dabij)["cdij"] += ( + 0.5 ) * Vijkl["ijij"];
  (*Dabij)["ccij"] += ( + 1.0 ) * Viajb["icic"];
  (*Dabij)["cdij"] += ( - 1.0 ) * Viajb["icic"];
  (*Dabij)["ccii"] += ( - 1.0 ) * Viajb["icic"];
  (*Dabij)["cdii"] += ( + 1.0 ) * Viajb["icic"];
  (*Dabij)["cdij"] += ( + 0.5 ) * Vabcd["cdcd"];
  (*Dabij)["ccij"] += ( + 0.5 ) * Tabij["ecij"] * Vijab["ijec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["ecij"] * Vijab["ijec"];
  (*Dabij)["cdij"] += ( + 0.25) * Tabij["efij"] * Vijab["ijef"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dabij)["cdii"] += ( + 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dabij)["ccij"] += ( - 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["ccii"] += ( + 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdii"] += ( - 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["efoi"] * Vijab["oief"];
  (*Dabij)["cdii"] += ( + 0.5 ) * Tabij["efoi"] * Vijab["oief"];
  (*Dabij)["cdij"] += ( + 0.25) * Tabij["cdmn"] * Vijab["mncd"];
  (*Dabij)["ccij"] += ( + 0.5 ) * Tabij["ecno"] * Vijab["noec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["ecno"] * Vijab["noec"];
*/
}

template <typename F>
std::vector<SDFockVector<F>>
CcsdPreconditioner<F>::getInitialBasis(const int eigenVectorsCount) {
  LOG(0, "CcsdPreconditioner") << "Getting initial basis " << std::endl;
  DefaultRandomEngine randomEngine;
  std::normal_distribution<double> normalDistribution(
    0.0, preconditionerRandomSigma
  );
  if (preconditionerRandom) {
    LOG(0, "CcsdPreconditioner") << "Randomizing initial guess" << std::endl;
  }
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements( diagonalH.readLocal() );
  std::sort(
    localElements.begin(), localElements.end(),
    EomDiagonalValueComparator<F>()
  );
  int localElementsSize( localElements.size() );

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10*eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  for (
    int i(0);
    i < std::min(localElementsSize, trialEigenVectorsCount);
    ++i
  ) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Cc4s::world);
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
  std::sort(
    lowestElements.begin(), lowestElements.end(),
    EomDiagonalValueComparator<F>()
  );
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;

  int currentEigenVectorCount(0);
  unsigned int b(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    V basisElement(diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t,F>> elements;
    if (communicator.getRank() == 0) {
      if ( b >= lowestElements.size() ) {
        throw EXCEPTION("No more elements to create initial basis");
      }
      elements.push_back(
        std::make_pair(lowestElements[b].first, 1.0)
      );
    }
    basisElement.write(elements);
    if (preconditionerRandom) {
      auto Rai(*basisElement.get(0));
      auto Rabij(*basisElement.get(1));
      setRandomTensor(Rai, normalDistribution, randomEngine);
      setRandomTensor(Rabij, normalDistribution, randomEngine);
      (*basisElement.get(0))["ai"] += Rai["ai"];
      (*basisElement.get(1))["abij"] += Rabij["abij"];
    }
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    (*basisElement.get(1))["abii"]=0.0;
    (*basisElement.get(1))["aaij"]=0.0;
    (*basisElement.get(1))["aaii"]=0.0;

    double preDot2(std::abs(basisElement.dot(basisElement)));
    // Antisymmetrize the new basis element
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["abji"];
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["baij"];

    OUT() << "\tnormPreSymmetrize=" << preDot2 << std::endl;

    double preDot3(std::abs(basisElement.dot(basisElement)));
    OUT() << "\tnormAfterSymmetrize=" << preDot3 << std::endl;

    OUT() << "\tbasisSize=" << basis.size() << std::endl;

    // Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basis[j].dot(basisElement);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if ( std::abs(basisElementNorm) < 1e-10 ) {
      zeroVectorCount++;
      b++;
      continue;
    }

    basisElement = 1.0 / std::sqrt(basisElement.dot(basisElement))*basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    b++;

    if ( std::abs(basisElementNorm - double(1)) > 1e-10 * double(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);

  }

  // Now make sure that you are giving an antisymmetric basis
  // and also that it is grammschmited afterwards
  //LOG(1,"CcsdPreconditioner") << "Antisymmetrize basis" << std::endl;
  //for (unsigned int j(0); j < basis.size(); ++j) {
    //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["abji"];
    //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["baij"];
  //}

  //LOG(1,"CcsdPreconditioner") <<
      //"Performing Gramm Schmidt in the initial basis"
  //<< std::endl;
  //for (unsigned int b(0); b < basis.size(); ++b) {
    //V newVector(basis[b]);
    //for (unsigned int j(0); j < b; ++j) {
      //newVector -= basis[j] * basis[j].dot(basis[b]);
    //}
    //// normalize
    //basis[b] = 1 / std::sqrt(newVector.dot(newVector)) * newVector;
  //}

    OUT() << "\treturning vecs =" << basis.size() << std::endl;

  return basis;
}

template <typename F>
CcsdtFockVector<F>
CcsdPreconditioner<F>::getCorrection(
  const complex lambda, CcsdtFockVector<F> &residuum
) {
  // Cast ccsdt into ccsd
  V w(residuum);
  // apply the old getCorrection
  // and cast again into a ccsdt vector
  CcsdtFockVector<F> result(getCorrection(lambda, w));
  // et voila
  return result;

}

template <typename F>
SDFockVector<F>
CcsdPreconditioner<F>::getCorrection(
  const complex lambda, SDFockVector<F> &residuum
) {
  SDFockVector<F> w(diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
    public:
      DiagonalCorrection(const double lambda_): lambda(lambda_) {
      }
      F operator ()(const F residuumElement, const F diagonalElement) {
        return std::abs(lambda - diagonalElement) < 1E-4 ?
          0.0 : residuumElement / (lambda - diagonalElement);
      }
    protected:
      double lambda;
  } diagonalCorrection(std::real(lambda));

  SDFockVector<F> correction(diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.getComponentsCount(); ++c) {
    const char *indices( correction.componentIndices[c].c_str() );
    (*correction.get(c)).contract(
      1.0,
      *residuum.get(c),indices,
      *diagonalH.get(c),indices,
      0.0,indices,
      CTF::Bivar_Function<F>(diagonalCorrection)
    );
  }
  // Filter out unphysical components from the correction
  (*correction.get(1))["abii"]=0.0;
  (*correction.get(1))["aaij"]=0.0;
  (*correction.get(1))["aaii"]=0.0;

  // Antisymmetrize the correction
  (*correction.get(1))["abij"] -= (*correction.get(1))["abji"];
  (*correction.get(1))["abij"] -= (*correction.get(1))["baij"];
  (*correction.get(1))["abij"] = 0.25 * (*correction.get(1))["abij"];

  return correction;
}

// instantiate
template class CcsdPreconditioner<double>;
template class CcsdPreconditioner<complex>;


//template <> class SDFockVector::CcsdPreconditioner<double>{};
//template
//CcsdPreconditioner<double>::CcsdPreconditioner();
