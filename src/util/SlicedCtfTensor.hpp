#ifndef SLICED_CTF_TENSOR_DEFINED
#define SLICED_CTF_TENSOR_DEFINED

#include <math/Permutation.hpp>
#include <util/Tensor.hpp>
#include <initializer_list>
#include <vector>
#include <string>

namespace sisi4s {
template <typename F = double>
class SlicedCtfTensor {
public:
  SlicedCtfTensor(Tensor<F> &T, const std::initializer_list<int> &list) {
    create(T, std::vector<int>(list));
  }
  template <typename Iterator>
  SlicedCtfTensor(Tensor<F> &T, Iterator begin, Iterator end) {
    create(T, std::vector<int>(begin, end));
  }
  void create(Tensor<F> &T, const std::vector<int> &slicedDims) {
    // generate indexNames = "abcd..."
    std::string indexNames(T.order, 'a');
    for (int d(0); d < T.order; ++d) indexNames[d] += d;
    // figure out how many slices we have in each sliced dimension
    slicedLens.resize(slicedDims.size());
    int64_t slicesCount = 1;
    for (unsigned int d(0); d < slicedDims.size(); ++d) {
      slicesCount *= slicedLens[d] = T.lens[slicedDims[d]];
    }
    // do all slices
    slices.resize(slicesCount);
    for (int64_t sliceIndex(0); sliceIndex < slicesCount; ++sliceIndex) {
      // get slice position from its index number
      int64_t i(sliceIndex);
      std::vector<int> slicePosition(slicedDims.size());
      for (unsigned int d(0); d < slicedDims.size(); ++d) {
        slicePosition[d] = i % slicedLens[d];
        i /= slicedLens[d];
      }
      // start = (0,...,0), end = (lens[0],...,lend[order-1])
      std::vector<int> start(T.order), end(T.lens, T.lens + T.order);
      // now enter the slice position in the sliced dimensions
      for (unsigned int d(0); d < slicedDims.size(); ++d) {
        start[slicedDims[d]] = slicePosition[d];
        end[slicedDims[d]] = slicePosition[d] + 1;
      }
      slices[sliceIndex] = new Tensor<F>(T.slice(start.data(), end.data()));
    }
  }
  ~SlicedCtfTensor() {
    for (uint64_t i(0); i < slices.size(); ++i) {
      if (slices[i]) delete slices[i];
      slices[i] = nullptr;
    }
  }

  Tensor<F> &operator()(const std::initializer_list<int> &slicePos) {
    return (*this)(std::vector<int>(slicePos));
  }
  template <typename Iterator>
  Tensor<F> &operator()(Iterator begin, Iterator end) {
    return (*this)(std::vector<int>(begin, end));
  }
  Tensor<F> &operator()(const std::vector<int> &slicePosition) {
    // get slice number from its position
    int64_t sliceIndex(0);
    for (int d(slicedLens.size() - 1); d >= 0; --d) {
      sliceIndex *= slicedLens[d];
      sliceIndex += slicePosition[d];
    }
    return *slices[sliceIndex];
  }

  std::vector<int> slicedLens;
  std::vector<Tensor<F> *> slices;
};
} // namespace sisi4s

#endif
