/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SRC_ALGORITHMS_COMPLEXTENSORREADER_HPP_
#define SRC_ALGORITHMS_COMPLEXTENSORREADER_HPP_
#include <vector>

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

class ComplexTensorReader : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(ComplexTensorReader);
  explicit ComplexTensorReader(std::vector<Argument> const &argumentList);
  virtual ~ComplexTensorReader();
  /**
   * \brief Reads real tensor data into the tensor Data.
   */
  virtual void run();
};

} // namespace sisi4s

#endif //  SRC_ALGORITHMS_COMPLEXTENSORREADER_HPP_
