/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/Algorithm.hpp>
#include <Data.hpp>
#include <math/Complex.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>

using namespace sisi4s;

Algorithm::Algorithm(std::vector<Argument> const &argumentList) {
  for (auto arg(argumentList.begin()); arg != argumentList.end(); ++arg) {
    Argument argument = *arg;
    arguments[argument.getName()] = argument.getData();
  }
}

Algorithm::~Algorithm() {
}

/**
 * \brief The dryRun estimates resource consumption, especially
 * memory and processor time.
 */
void Algorithm::dryRun() {
  LOG(0, getName()) << "dry run not implemented" << std::endl;
}

bool Algorithm::isArgumentGiven(std::string const &name) {
  return arguments.find(name) != arguments.end();
}

Data *Algorithm::getArgumentData(std::string const &name) {
  auto dataIterator(arguments.find(name));
  if (dataIterator == arguments.end()) {
    std::stringstream sStream;
    sStream << "Missing argument: " << name;
//    throw new EXCEPTION(std::stringstream() << "Missing argument: " << name);
    throw new EXCEPTION(sStream.str());
  }
  Data *data = Data::get(dataIterator->second);
  if (!data) {
    std::stringstream sStream;
    sStream << "Missing data: " << dataIterator->second;
//    throw new EXCEPTION(std::stringstream() << "Missing data: " << dataIterator->second);
    throw new EXCEPTION(sStream.str());
  }
  return data;
}

std::string Algorithm::getTextArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  TextData const *textData = dynamic_cast<TextData const *>(data);
  if (!textData) {
    std::stringstream sstream;
    sstream << "Incompatible type for argument: " << name << ". "
      << "Excpected Text, found " << data->getTypeName() << ".";
    throw new EXCEPTION(sstream.str());
  }
  return textData->value;
}
std::string Algorithm::getTextArgument(
  std::string const &name, std::string const &defaultValue
) {
  return isArgumentGiven(name) ? getTextArgument(name) : defaultValue;
}

bool Algorithm::getBooleanArgument(std::string const &name) {
  /*
   *TODO: Do this without the getTextArgument function, because in this
   *case the parser want to have quotes on the boolean value i.e.
   *  (myflag "true")
   */
  std::string text(getTextArgument(name));
  if (
    text.compare(".TRUE.") == 0 ||
    text.compare("true") == 0 ||
    text.compare("True") == 0 ||
    text.compare("TRUE") == 0 ||
    text.compare("1") == 0 ||
    text.compare("t") == 0 ||
    text.compare("T") == 0
  ) {
    return true;
  } else {
    return false;
  }
}
bool Algorithm::getBooleanArgument(
  std::string const &name, bool const &defaultValue
) {
  return isArgumentGiven(name) ? getBooleanArgument(name) : defaultValue;
}

int64_t Algorithm::getIntegerArgument(std::string const &name) {
  Data const *data(getArgumentData(name));
  IntegerData const *integerData = dynamic_cast<IntegerData const *>(data);
  if (!integerData) {
    std::stringstream sstream;
    sstream << "Incompatible type for argument: " << name << ". "
      << "Excpected Integer, found " << data->getTypeName() << ".";
    throw new EXCEPTION(sstream.str());
  }
  return integerData->value;
}
int64_t Algorithm::getIntegerArgument(
  std::string const &name, int64_t const defaultValue
) {
  return isArgumentGiven(name) ? getIntegerArgument(name) : defaultValue;
}

sisi4s::real Algorithm::getRealArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  RealData *realData(dynamic_cast<RealData *>(data));
  if (realData) return realData->value;
  IntegerData *integerData(dynamic_cast<IntegerData *>(data));
  if (integerData) return getRealArgumentFromInteger(integerData);
  TensorData<real> *tensorData(dynamic_cast<TensorData<real> *>(data));
  if (tensorData) return getRealArgumentFromTensor(tensorData);
  std::stringstream sstream;
  sstream << "Incompatible type for argument: " << name << ". "
    << "Excpected Real, found " << data->getTypeName() << ".";
  throw new EXCEPTION(sstream.str());
}
sisi4s::real Algorithm::getRealArgument(
  const std::string &name, const real defaultValue
) {
  return isArgumentGiven(name) ? getRealArgument(name) : defaultValue;
}
sisi4s::real Algorithm::getRealArgumentFromInteger(IntegerData *integerData ) {
  real value(integerData->value);
  if (int64_t(value) != integerData->value) {
    LOG(0, "root") << "Warning: loss of precision in conversion from integer to real."
      << std::endl;
  }
  return value;
}
sisi4s::real Algorithm::getRealArgumentFromTensor(TensorData<real> *data) {
  Assert(
    data->value->order == 0,
    "Scalar expected in conversion from tensor to real."
  );
  // retrieve the real value from the tensor
  CTF::Scalar<real> scalar;
  scalar[""] = (*data->value)[""];
  return scalar.get_val();
}

template <typename F, typename C>
C *Algorithm::getContainerArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  ContainerData<F, C> *containerData(dynamic_cast<ContainerData<F, C> *>(data));
  if (containerData) return containerData->value;
  // TODO: provide conversion routines from real to complex containers
  std::stringstream sStream;
  sStream << "Incompatible type for argument: " << name << ". "
    << "Excpected container of " << TypeTraits<F>::getName()
    << ", found " << data->getTypeName() << ".";
  throw new EXCEPTION(sStream.str());
}

template <typename F, typename C>
void Algorithm::allocateContainerArgument(
  std::string const &name, C *container
) {
  Data *mentionedData(getArgumentData(name));
  new ContainerData<F, C >(
    mentionedData->getName(), container
  );
}
template
void Algorithm::allocateContainerArgument< int64_t, std::vector<int64_t> >(
  std::string const &name, std::vector<int64_t> *container);
template
void Algorithm::allocateContainerArgument< sisi4s::real, std::vector<sisi4s::real> >(
  std::string const &name, std::vector<sisi4s::real> *container);
template
void Algorithm::allocateContainerArgument< sisi4s::complex, std::vector<sisi4s::complex> >(
  std::string const &name, std::vector<sisi4s::complex> *container);


template <typename F, typename T>
T *Algorithm::getTensorArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  TensorData<F, T> *tensorData(dynamic_cast<TensorData<F, T> *>(data));
  if (tensorData) return tensorData->value;
  RealData *realData(dynamic_cast<RealData *>(data));
  if (realData) return getTensorArgumentFromReal<F, T>(realData);
  // TODO: provide conversion routines from real to complex tensors
  std::stringstream sStream;
  sStream << "Incompatible type for argument: " << name << ". "
    << "Excpected tensor of " << TypeTraits<F>::getName()
    << ", found " << data->getTypeName() << ".";
  throw new EXCEPTION(sStream.str());
}
// instantiate
template
Tensor<Float64> *Algorithm::getTensorArgument<
  Float64, Tensor<Float64>
>(std::string const &);
template
Tensor<Complex64> *Algorithm::getTensorArgument<
  Complex64, Tensor<Complex64>
>(std::string const &);
template
DryTensor<Float64> *Algorithm::getTensorArgument<
  Float64, DryTensor<Float64>
>(std::string const &);
template
DryTensor<Complex64> *Algorithm::getTensorArgument<
  Complex64, DryTensor<Complex64>
>(std::string const &);


/**
 * \brief Traits for retrieving the Scalar, Vector and Matrix tensor type.
 */
template < typename F, typename T=Tensor<F> >
class TensorTypeTraits;

template <typename F>
class TensorTypeTraits< F, Tensor<F> > {
public:
  typedef Tensor<F> BaseType;
  typedef CTF::Scalar<F> ScalarType;
  typedef CTF::Vector<F> VectorType;
  typedef CTF::Matrix<F> MatrixType;
};
template <typename F>
class TensorTypeTraits< F, CTF::Matrix<F> > {
public:
  typedef Tensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, CTF::Vector<F> > {
public:
  typedef Tensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, CTF::Scalar<F> > {
public:
  typedef Tensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, DryTensor<F> > {
public:
  typedef DryTensor<F> BaseType;
  typedef DryScalar<F> ScalarType;
  typedef DryVector<F> VectorType;
  typedef DryMatrix<F> MatrixType;
};
template <typename F>
class TensorTypeTraits< F, DryMatrix<F> > {
public:
  typedef DryTensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, DryVector<F> > {
public:
  typedef DryTensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, DryScalar<F> > {
public:
  typedef DryTensor<F> BaseType;
};


/**
 * \brief Converts the given real data into a scalar tensor.
 */
template <typename F, typename T>
T *Algorithm::getTensorArgumentFromReal(RealData *realData) {
  // FIXME: left to leak memory...
  // a better solution would be to replace the RealData with the allocated
  // TensorData and support down-cast for Scalars to Real
  return new typename TensorTypeTraits<F,T>::ScalarType(realData->value);
}
// instantiate
template
Tensor<Float64> *Algorithm::getTensorArgumentFromReal<
  Float64, Tensor<Float64>
>(RealData *);
template
Tensor<Complex64> *Algorithm::getTensorArgumentFromReal<
  Complex64, Tensor<Complex64>
>(RealData *);
template
DryTensor<Float64> *Algorithm::getTensorArgumentFromReal<
  Float64, DryTensor<Float64>
>(RealData *);
template
DryTensor<Complex64> *Algorithm::getTensorArgumentFromReal<
  Complex64, DryTensor<Complex64>
>(RealData *);

template <typename F, typename T>
void Algorithm::allocatedTensorArgument(
  std::string const &name, T *tensor
) {
  Data *mentionedData(getArgumentData(name));
  new TensorData<F, typename TensorTypeTraits<F, T>::BaseType>(
    mentionedData->getName(), tensor
  );
  // NOTE: the constructor of TensorData enteres its location in the
  // data map and destroys the previous content, i.e. mentionedData.
}
// instantiate
template
void Algorithm::allocatedTensorArgument<
  Float64, Tensor<Float64>
>(std::string const &name, Tensor<Float64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
template
void Algorithm::allocatedTensorArgument<
  sisi4s::real, CTF::Matrix<sisi4s::real>
>(std::string const &name, CTF::Matrix<real> *tensor);
template
void Algorithm::allocatedTensorArgument<
  sisi4s::real, CTF::Vector<sisi4s::real>
>(std::string const &name, CTF::Vector<sisi4s::real> *tensor);
template
void Algorithm::allocatedTensorArgument<
  sisi4s::real, CTF::Scalar<sisi4s::real>
>(std::string const &name, CTF::Scalar<sisi4s::real> *tensor);

template
void Algorithm::allocatedTensorArgument<
  Complex64, Tensor<Complex64>
>(std::string const &name, Tensor<Complex64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
template
void Algorithm::allocatedTensorArgument<
  sisi4s::complex, CTF::Matrix<sisi4s::complex>
>(std::string const &name, CTF::Matrix<sisi4s::complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  sisi4s::complex, CTF::Vector<sisi4s::complex>
>(std::string const &name, CTF::Vector<sisi4s::complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  sisi4s::complex, CTF::Scalar<sisi4s::complex>
>(std::string const &name, CTF::Scalar<sisi4s::complex> *tensor);

template
void Algorithm::allocatedTensorArgument<
  Float64, DryTensor<Float64>
>(std::string const &name, DryTensor<Float64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
template
void Algorithm::allocatedTensorArgument<
  Float64, DryMatrix<Float64>
>(std::string const &name, DryMatrix<Float64> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Float64, DryVector<Float64>
>(std::string const &name, DryVector<Float64> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Float64, DryScalar<Float64>
>(std::string const &name, DryScalar<Float64> *tensor);

template
void Algorithm::allocatedTensorArgument<
  Complex64, DryTensor<Complex64>
>(std::string const &name, DryTensor<Complex64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
template
void Algorithm::allocatedTensorArgument<
  Complex64, DryMatrix<Complex64>
>(std::string const &name, DryMatrix<Complex64> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Complex64, DryVector<Complex64>
>(std::string const &name, DryVector<Complex64> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Complex64, DryScalar<Complex64>
>(std::string const &name, DryScalar<Complex64> *tensor);


void Algorithm::setRealArgument(std::string const &name, const real value) {
  Data *mentionedData(getArgumentData(name));
  new RealData(mentionedData->getName(), value);
}

void Algorithm::setIntegerArgument(std::string const &name, int const value) {
  Data *mentionedData(getArgumentData(name));
  new IntegerData(mentionedData->getName(), value);
}

AlgorithmFactory::AlgorithmMap *AlgorithmFactory::algorithmMap;

