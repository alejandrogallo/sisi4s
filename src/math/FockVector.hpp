#ifndef FOCK_VECTOR_DEFINED
#define FOCK_VECTOR_DEFINED

#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <math/MathFunctions.hpp>
#include <Sisi4s.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <ostream>

#include <util/Tensor.hpp>

namespace sisi4s {
template <typename F = double>
/**
 * \brief Represents the direct sum of Tensors and provides the
 * vector space operations of addition, scalar multiplication, inner product,
 * complex conjugation to get dual vectors and matrix multiplication
 * between vectors and duals, which yields a scalar.
 **/
class FockVector {
public:
  typedef F FieldType;

  std::vector<PTR(Tensor<F>)> component_tensors;
  std::vector<std::string> component_indices;

  /**
   * \brief Default constructor for an empty Fock vector without elements.
   **/
  FockVector() {}

  /**
   * \brief Move constructor taking possession of the tensors owned by a.
   **/
  FockVector(FockVector<F> &&a)
      : component_tensors(a.component_tensors)
      , component_indices(a.component_indices)
      , index_ends(a.component_tensors.size()) {
    build_index_translation();
  }

  /**
   * \brief Copy constructor copying the tensors owned by a.
   **/
  FockVector(const FockVector<F> &a)
      : component_tensors(a.component_tensors.size())
      , component_indices(a.component_indices)
      , index_ends(a.component_tensors.size()) {
    copy_components(a.component_tensors);
    build_index_translation();
  }

  /**
   * \brief Move constructor taking possession of the tensors given.
   **/
  FockVector(const std::vector<PTR(Tensor<F>)> &tensors,
             const std::vector<std::string> &indices)
      : component_tensors(tensors)
      , component_indices(indices)
      , index_ends(component_tensors.size()) {
    build_index_translation();
  }

  /**
   * \brief Move constructor taking possession of the tensors given
   * by the iterators.
   **/
  template <typename TensorsIterator, typename IndicesIterator>
  FockVector(TensorsIterator tensorsBegin,
             TensorsIterator tensorsEnd,
             IndicesIterator indicesBegin,
             IndicesIterator indicesEnd)
      : component_tensors(tensorsBegin, tensorsEnd)
      , component_indices(indicesBegin, indicesEnd)
      , index_ends(component_tensors.size()) {
    build_index_translation();
  }

  /**
   * \brief Retrieves the i-th component tensor. Note that
   * the Tensor is not const since rearrangement may be
   * required also in non-modifying tensor operations.
   **/
  const PTR(Tensor<F>) &get(const size_t i) const {
    return component_tensors[i];
  }

  /**
   * \brief Retrieves the i-th component tensor.
   **/
  PTR(Tensor<F>) &get(const size_t i) { return component_tensors[i]; }

  /**
   * \brief Retrieves the i-th component indices.
   **/
  const std::string &get_indices(const size_t i) const {
    return component_indices[i];
  }

  /**
   * \brief Retrieves the i-th component indices as modifiable string.
   **/
  std::string &get_indices(const size_t i) { return component_indices[i]; }

  /**
   * \brief Move assignment operator taking possession of the tensors
   * owned by a.
   **/
  FockVector<F> &operator=(const FockVector<F> &&a) {
    component_tensors = a.component_tensors;
    component_indices = a.component_indices;
    build_index_translation();
    return *this;
  }

  /**
   * \brief Copy assignment operator copying the tensors owned by a.
   **/
  FockVector<F> &operator=(const FockVector<F> &a) {
    component_indices = a.component_indices;
    copy_components(a.component_tensors);
    build_index_translation();
    return *this;
  }

  /**
   * \brief Add-to assignment operator adding each component of a
   * to the respective component of this FockVector.
   **/
  FockVector<F> &operator+=(const FockVector<F> &a) {
    checkCompatibilityTo(a);
    for (size_t i(0); i < component_tensors.size(); ++i) {
      const char *indices(component_indices[i].c_str());
      get(i)->sum(+1.0, *a.get(i), indices, 1.0, indices);
    }
    return *this;
  }

  /**
   * \brief Subtract-from assignment operator subtracting each component of a
   * from the respective component of this FockVector.
   **/
  FockVector<F> &operator-=(const FockVector<F> &a) {
    checkCompatibilityTo(a);
    for (size_t i(0); i < component_tensors.size(); ++i) {
      const char *indices(get_indices(i).c_str());
      get(i)->sum(-1.0, *a.get(i), indices, 1.0, indices);
    }
    return *this;
  }

  /**
   * \brief Multiply-by assignment operator scalar multiplying each component
   * each component of this FockVector by the given scalar.
   **/
  FockVector<F> &operator*=(const F s) {
    for (size_t i(0); i < component_tensors.size(); ++i) {
      const char *indices(get_indices(i).c_str());
      get(i)->sum(s, *get(i), indices, 0.0, indices);
    }
    return *this;
  }

  /**
   * \brief Creates and returns the conjugate transpose of this FockVector.
   * The first and the second half of the inidices in each component are
   * swapped for the transposition. For real types F the conjugation
   * does nothing.
   **/
  FockVector<F> conjugateTranspose() const {
    FockVector<F> result;
    for (size_t i(0); i < component_tensors.size(); ++i) {
      size_t order(get_indices(i).length() / 2);
      std::vector<int> transposedLens(get(i)->lens, get(i)->lens + 2 * order);
      std::rotate(transposedLens.begin(),
                  transposedLens.begin() + order,
                  transposedLens.begin() + 2 * order);
      result.component_tensors.push_back(
          NEW(Tensor<F>,
              transposedLens.size(),
              transposedLens.data(),
              get(i)->sym,
              *get(i)->wrld,
              (std::string(get(i)->get_name()) + "*").c_str()));
      result.component_indices.push_back(get_indices(i).substr(order, 2 * order)
                                         + get_indices(i).substr(0, order));
      CTF::Univar_Function<F> fConj(sisi4s::conj<F>);
      result.get(i)->sum(1.0,
                         *get(i),
                         get_indices(i).c_str(),
                         0.0,
                         result.get_indices(i).c_str(),
                         fConj);
    }
    return result;
  }

  /**
   * \brief Returns the matrix product of this bra-FockVector with the
   * given dual ket-FockVector ket.
   **/
  F braket(const FockVector<F> &ket) const {
    // checkDualCompatibility(ket);
    CTF::Scalar<F> result;
    for (size_t i(0); i < component_tensors.size(); ++i) {
      const char *indices(get_indices(i).c_str());
      const char *ketIndices(ket.get_indices(i).c_str());
      // add to result
      result[""] += (*get(i))[indices] * (*ket.get(i))[ketIndices];
    }
    return result.get_val();
  }

  /**
   * \brief Returns the inner product of this ket-FockVector with the
   * given ket-FockVector a. The elements of this FockVector are conjugated
   * in the inner product, i.e. this->dot(a) yields the same results as
   * this->conjugateTranspose().braket(a).
   **/
  F dot(const FockVector<F> &a) const {
    checkCompatibilityTo(a);
    CTF::Scalar<F> result;
    for (size_t i(0); i < component_tensors.size(); ++i) {
      const char *indices(get_indices(i).c_str());
      CTF::Bivar_Function<F> fDot(&sisi4s::dot<F>);
      // add to result
      result.contract(1.0, *get(i), indices, *a.get(i), indices, 1.0, "", fDot);
    }
    return result.get_val();
  }

  /**
   * \brief Get the number of component tensors of this FockVector.
   */
  size_t get_components_count() const { return component_tensors.size(); }

  /**
   * \brief Get the total number of degrees of freedom represented by this
   * FockVector, i.e. the total number of field values contained in all
   * component tensors. The indices used by read and write are between
   * 0 and get_dimension()-1.
   */
  size_t get_dimension() const { return index_ends.back(); }

  /**
   * \Brief Translates the given component and component index into
   * its element into an index between 0 and get_dimension()-1.
   **/
  size_t get_index(const size_t component, const size_t componentIndex) const {
    size_t base(component > 0 ? index_ends[component - 1] : 0);
    return base + componentIndex;
  }

  /**
   * \Brief Translates the given index between 0 and get_dimension()-1
   * into a component number and component index into the corresponding
   * component tensor.
   **/
  void fromIndex(const size_t index,
                 size_t &component,
                 size_t &componentIndex) const {
    component = 0;
    size_t base(0);
    while (component < index_ends.size()) {
      if (index < index_ends[component]) break;
      base = index_ends[component++];
    }
    if (component >= index_ends.size()) {
      throw new EXCEPTION("Index out bounds");
    }
    componentIndex = index - base;
  }

  /**
   * \brief Reads out all locally stored values together with their
   * respective indices. The indices are between 0 and get_dimension()-1.
   **/
  std::vector<std::pair<size_t, F>> readLocal() const {
    size_t elementsCount(0);
    std::vector<std::pair<size_t, F>> elements;
    for (size_t i(0); i < component_tensors.size(); ++i) {
      size_t componentValuesCount;
      size_t *component_indices;
      F *componentValues;
      get(i)->read_local(reinterpret_cast<int64_t *>(&componentValuesCount),
                         reinterpret_cast<int64_t **>(&component_indices),
                         &componentValues);

      elements.resize(elementsCount + componentValuesCount);
      for (size_t k(0); k < componentValuesCount; ++k) {
        // translate index within component tensor to FockVector index
        elements[elementsCount + k].first = get_index(i, component_indices[k]);
        elements[elementsCount + k].second = componentValues[k];
      }
      elementsCount += componentValuesCount;
      free(component_indices);
      free(componentValues);
    }
    return elements;
  }

  /**
   * \brief Writes the given values together with their
   * respective indices. The indices are between 0 and get_dimension()-1.
   **/
  void write(const std::vector<std::pair<size_t, F>> &elements) {
    // vectors to contain indices and values for each component tensor
    std::vector<std::vector<size_t>> tensorIndices(component_tensors.size());
    std::vector<std::vector<F>> tensorValues(component_tensors.size());

    for (size_t k(0); k < elements.size(); ++k) {
      size_t component;
      size_t componentIndex;
      fromIndex(elements[k].first, component, componentIndex);
      // write to respective component tensor
      tensorIndices[component].push_back(componentIndex);
      tensorValues[component].push_back(elements[k].second);
    }

    // write data of each tensor
    for (size_t i(0); i < component_tensors.size(); ++i) {
      tensorIndices[i].reserve(tensorIndices[i].size() + 1);
      tensorValues[i].reserve(tensorIndices[i].size() + 1);
      get(i)->write(tensorIndices[i].size(),
                    reinterpret_cast<int64_t *>(tensorIndices[i].data()),
                    tensorValues[i].data());
    }
  }

protected:
  /**
   * \brief The end of the FockVector index range for each component.
   * This vector is used for translating component number and indices
   * into FockVector indicies.
   **/
  std::vector<size_t> index_ends;

  /**
   * \Brief Builds the index ends vector needed for the
   * index translation methods get_index and fromIndex.
   **/
  void build_index_translation() {
    index_ends.resize(component_tensors.size());
    size_t indexBase(0);
    for (size_t i(0); i < component_tensors.size(); ++i) {
      size_t tensorIndexSize(1);
      for (int d(0); d < get(i)->order; ++d) {
        tensorIndexSize *= get(i)->lens[d];
      }
      index_ends[i] = indexBase += tensorIndexSize;
    }
  }

  /**
   * \brief Sets this FockVector's component tensors by copying the given
   * component tensors. Called by copy constructors and copy assignments.
   **/
  void copy_components(const std::vector<PTR(Tensor<F>)> &components) {
    component_tensors.resize(components.size());
    for (size_t i(0); i < components.size(); ++i) {
      component_tensors[i] = NEW(Tensor<F>, *components[i]);
    }
  }

  /**
   * \brief Check if two FockVectors are transpose of each other by swapping
   * the first and the second half of the component indices.
   **/
  // TODO: Improve speed?
  void checkDualCompatibility(const FockVector<F> &a) const {
    checkCompatibilityTo(a);
    for (size_t i(0); i < component_tensors.size(); i++) {
      size_t indexLens(a.get(i)->order());
      for (size_t j(0); j < indexLens; j++) {
        size_t indexPos(get(i).find(a.getIndicies(i)[j]));
        if (indexPos == std::string::npos) {
          throw EXCEPTION("Indices of fock vectors do not match");
        }
        if (a.get(i)->lens[j] != get(i)->lens[indexPos]) {
          throw EXCEPTION("Shapes of component tensors does not match");
        }
      }
    }
  }

  void checkCompatibilityTo(const FockVector<F> &a) const {
    if (component_tensors.size() != a.component_tensors.size()
        || component_indices.size() != a.component_indices.size()) {
      throw EXCEPTION("Number of component tensors does no match");
    }
    // TODO: check shapes.
  }
};

/**
 * \brief Returns the sum of two FockVectors a and b, where
 * neither a nor b are modified.
 **/
template <typename F>
inline FockVector<F> operator+(const FockVector<F> &a, const FockVector<F> &b) {
  FockVector<F> result(a);
  result += b;
  return result;
}
/**
 * \brief Returns the sum of two FockVectors a and b, where
 * a is movable and will be used for the result.
 **/
template <typename F>
inline FockVector<F> &&operator+(FockVector<F> &&a, const FockVector<F> &b) {
  a += b;
  return a;
}
/**
 * \brief Returns the sum of two FockVectors a and b, where
 * b is movable and will be used for the result.
 **/
template <typename F>
inline FockVector<F> &&operator+(FockVector<F> &a, const FockVector<F> &&b) {
  b += a;
  return b;
}

/**
 * \brief Returns the difference between two FockVectors a and b, where
 * neither a nor b are modified.
 **/
template <typename F>
inline FockVector<F> operator-(const FockVector<F> &a, const FockVector<F> &b) {
  FockVector<F> result(a);
  result -= b;
  return result;
}
/**
 * \brief Returns the difference between two FockVectors a and b, where
 * a is movable and will be used for the result.
 **/
template <typename F>
inline FockVector<F> &&operator-(FockVector<F> &&a, const FockVector<F> &b) {
  a -= b;
  return a;
}
/**
 * \brief Returns the difference between two FockVectors a and b, where
 * b is movable and will be used for the result.
 **/
template <typename F>
inline FockVector<F> &&operator-(const FockVector<F> &a, FockVector<F> &&b) {
  b -= a;
  // TODO: directly invoke sum to prevent extra multiplication by -1
  b *= F(-1);
  return b;
}

/**
 * \brief Returns the scalar multiple of the FockVector a
 * right-multiplied with the scalar s, where a is not modified.
 **/
template <typename F>
inline FockVector<F> operator*(const FockVector<F> &a, const F s) {
  FockVector<F> result(a);
  result *= s;
  return result;
}
/**
 * \brief Returns the scalar multiple of the FockVector a
 * right-multiplied with the scalar s, where a movable and will be used
 * for the result.
 **/
template <typename F>
inline FockVector<F> &&operator*(FockVector<F> &&a, const F s) {
  a *= s;
  return a;
}

/**
 * \brief Returns the scalar multiple of the FockVector a
 * left-multiplied with the scalar s, where a is not modified.
 **/
template <typename F>
inline FockVector<F> operator*(const F s, const FockVector<F> &a) {
  FockVector<F> result(a);
  result *= s;
  return result;
}
/**
 * \brief Returns the scalar multiple of the FockVector a
 * left-multiplied with the scalar s, where a movable and will be used
 * for the result.
 **/
template <typename F>
inline FockVector<F> &&operator*(const F s, FockVector<F> &&a) {
  a *= s;
  return a;
}

/**
 * \brief Writes the FockVector a to the given stream and returns it
 * for further stream operations.
 **/
template <typename F>
inline std::ostream &operator<<(std::ostream &stream, const FockVector<F> &a) {
  stream << "( ";
  stream << a.get(0) << "[" << a.get_indices(0) << "]";
  for (size_t i(1); i < a.component_tensors.size(); ++i) {
    stream << ", " << a.get(i) << "[" << a.get_indices(i) << "]";
  }
  return stream << " )";
}

template <typename F, int N, int StartDimension = 0>
class FockVectorNdCanonical : public FockVector<F> {
public:
  using FockVector<F>::FockVector;

  /**
   * \brief Build a canonical vector from No and Nv.
   */
  FockVectorNdCanonical(const unsigned int No, const unsigned int Nv) {
    if (N > 6) {
      throw new EXCEPTION("FockVectorNdCanonical implemented only up to 6");
    }
    const std::string pindices("abcdefgABCDEFG");
    const std::string hindices("ijklomnIJKLOMN");
    for (unsigned int i(StartDimension); i <= N / 2; i++) {

      // vec<int> and not unsigned int, otherwise you'll be miserable for at
      // least two days.
      std::vector<int> syms(i, NS);
      std::vector<int> dimso(i, No);
      std::vector<int> dimsv(i, Nv);
      std::vector<int> dims(dimsv);
      dims.insert(dims.end(), dimso.begin(), dimso.end());

      this->component_indices.push_back(pindices.substr(0, i)
                                        + hindices.substr(0, i));
      this->component_tensors.push_back(
          NEW(Tensor<F>, 2 * i, dims.data(), syms.data(), *Sisi4s::world));
    }
    this->build_index_translation();
  }

  /**
   * \brief Empty constructor
   */
  FockVectorNdCanonical() {}

  /**
   * \brief Move constructor taking possession of the tensors owned by a.
   **/
  FockVectorNdCanonical(FockVector<F> &&a) {
    this->component_tensors = a.component_tensors;
    this->component_indices = a.component_indices;
    this->index_ends.resize(a.component_tensors.size());

    this->build_index_translation();
  }

  /**
   * \brief Copy constructor copying the tensors owned by a.
   **/
  FockVectorNdCanonical(const FockVector<F> &a) {
    this->component_tensors.resize(a.component_tensors.size());
    this->component_indices = a.component_indices;
    this->index_ends.resize(a.component_tensors.size());
    this->copy_components(a.component_tensors);

    this->build_index_translation();
  }
};

template <typename F>
using CISFockVector = FockVectorNdCanonical<F, 1, 0>;
template <typename F>
using CISDFockVector = FockVectorNdCanonical<F, 2, 0>;
template <typename F>
using CISDTFockVector = FockVectorNdCanonical<F, 3, 0>;

template <typename F>
class SDFockVector;

template <typename F>
class SFockVector : public FockVectorNdCanonical<F, 1, 1> {
public:
  using FockVectorNdCanonical<F, 1, 1>::FockVectorNdCanonical;

  SFockVector()
      : FockVectorNdCanonical<F, 1, 1>(0, 0) {}

  SFockVector(const SFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->index_ends.resize(1);
    this->component_tensors.resize(1);
    this->copy_components(a.component_tensors);
    this->build_index_translation();
  }

  SFockVector(SFockVector<F> &&a) {
    this->component_indices = a.component_indices;
    this->component_tensors = a.component_tensors;
    this->index_ends.resize(1);
    this->build_index_translation();
  }

  SFockVector<F> &operator=(const SFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->copy_components(a.component_tensors);
    this->build_index_translation();
    return *this;
  }

  SFockVector(const SDFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->index_ends.resize(1);
    this->component_tensors.resize(1);
    this->copy_components(a.component_tensors);
    this->build_index_translation();
  }
};

template <typename F>
class SDTFockVector;

template <typename F>
class SDFockVector : public FockVectorNdCanonical<F, 2, 1> {
public:
  using FockVectorNdCanonical<F, 2, 1>::FockVectorNdCanonical;

  SDFockVector()
      : FockVectorNdCanonical<F, 2, 1>(0, 0) {}

  SDFockVector(const SDFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->index_ends.resize(2);
    this->component_tensors.resize(2);
    this->copy_components(a.component_tensors);
    this->build_index_translation();
  }

  SDFockVector(SDFockVector<F> &&a) {
    this->component_indices = a.component_indices;
    this->component_tensors = a.component_tensors;
    this->index_ends.resize(2);
    this->build_index_translation();
  }

  SDFockVector<F> &operator=(const SDFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->copy_components(a.component_tensors);
    this->build_index_translation();
    return *this;
  }

  SDFockVector(const SFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->index_ends.resize(2);
    this->component_tensors.resize(2);
    this->copy_components(a.component_tensors);
    this->build_index_translation();
  }

  SDFockVector(const SDTFockVector<F> &a) {
    this->component_indices = a.component_indices;
    this->index_ends.resize(2);
    this->component_tensors.resize(2);
    this->copy_components(a.component_tensors);
    this->build_index_translation();
  }
};

template <typename F>
class SDTFockVector : public FockVectorNdCanonical<F, 3, 1> {
public:
  using FockVectorNdCanonical<F, 3, 1>::FockVectorNdCanonical;

  SDTFockVector()
      : FockVectorNdCanonical<F, 3, 1>(0, 0) {}

  SDTFockVector(const SDFockVector<F> &a) {
    this->copy_components(a.component_tensors);
    this->component_tensors.resize(3);
    this->component_indices.resize(3);
    this->component_indices[0] = a.component_indices[0];
    this->component_indices[1] = a.component_indices[1];
    this->component_indices[2] = "abcijk";
    this->index_ends.resize(3);
    // This copies the components of a

    int No(this->component_tensors[0]->lens[1]);
    int Nv(this->component_tensors[0]->lens[0]);
    int vvvooo[6] = {Nv, Nv, Nv, No, No, No};
    int syms[6] = {NS, NS, NS, NS, NS, NS};
    this->component_tensors[2] =
        NEW(Tensor<F>, 6, vvvooo, syms, *Sisi4s::world);
    (*this->get(2))["abcijk"] = 0.0;

    this->build_index_translation();
  }

  SDTFockVector(const SDFockVector<F> &&a) {
    this->component_tensors.resize(3);
    this->component_indices.resize(3);
    this->component_indices[0] = a.component_indices[0];
    this->component_indices[1] = a.component_indices[1];
    this->component_indices[2] = "abcijk";
    this->index_ends.resize(3);
    // This copies the components of a
    this->component_tensors[0] = a.component_tensors[0];
    this->component_tensors[1] = a.component_tensors[1];

    int No(this->component_tensors[0]->lens[1]);
    int Nv(this->component_tensors[0]->lens[0]);
    int vvvooo[6] = {Nv, Nv, Nv, No, No, No};
    int syms[6] = {NS, NS, NS, NS, NS, NS};
    this->component_tensors[2] =
        NEW(Tensor<F>, 6, vvvooo, syms, *Sisi4s::world);
    (*this->get(2))["abcijk"] = 0.0;

    this->build_index_translation();
  }

  SDTFockVector<F> &operator=(SDFockVector<F> &&a) {
    this->component_tensors.resize(3);
    this->component_indices.resize(3);
    this->component_indices[0] = a.component_indices[0];
    this->component_indices[1] = a.component_indices[1];
    this->component_indices[2] = "abcijk";
    this->index_ends.resize(3);
    // This copies the components of a
    this->component_tensors[0] = a.component_tensors[0];
    this->component_tensors[1] = a.component_tensors[1];

    int No(this->component_tensors[0]->lens[1]);
    int Nv(this->component_tensors[0]->lens[0]);
    int vvvooo[6] = {Nv, Nv, Nv, No, No, No};
    int syms[6] = {NS, NS, NS, NS, NS, NS};
    this->component_tensors[2] =
        NEW(Tensor<F>, 6, vvvooo, syms, *Sisi4s::world);
    (*this->get(2))["abcijk"] = 0.0;

    this->build_index_translation();

    return *this;
  }

  SDTFockVector<F> &operator=(const SDFockVector<F> &a) {
    this->copy_components(a.component_tensors);
    this->component_tensors.resize(3);
    this->component_indices.resize(3);
    this->component_indices[0] = a.component_indices[0];
    this->component_indices[1] = a.component_indices[1];
    this->component_indices[2] = "abcijk";
    this->index_ends.resize(3);
    // This copies the components of a

    int No(this->component_tensors[0]->lens[1]);
    int Nv(this->component_tensors[0]->lens[0]);
    int vvvooo[6] = {Nv, Nv, Nv, No, No, No};
    int syms[6] = {NS, NS, NS, NS, NS, NS};
    this->component_tensors[2] =
        NEW(Tensor<F>, 6, vvvooo, syms, *Sisi4s::world);
    (*this->get(2))["abcijk"] = 0.0;

    this->build_index_translation();

    return *this;
  }
};

} // namespace sisi4s

#endif
