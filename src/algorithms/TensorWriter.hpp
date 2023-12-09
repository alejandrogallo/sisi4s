#ifndef TENSOR_WRITER_DEFINED
#define TENSOR_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(TensorWriter,

                        template <typename F>
                        static void write(const std::string &name,
                                          const std::string fileName,
                                          Tensor<F> *A,
                                          const bool binary_p,
                                          const std::string rowIndexOrder,
                                          const std::string columnIndexOrder,
                                          const std::string delimiter);

);

} // namespace sisi4s

#endif
