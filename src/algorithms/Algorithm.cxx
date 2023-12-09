#include <iostream>
#include <sstream>
#include <vector>
#include <list>

#include <algorithms/Algorithm.hpp>
#include <NewData.hpp>
#include <Data.hpp>
#include <math/Complex.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>

using namespace sisi4s;

Algorithm::Algorithm(Arguments const &in_, Arguments const &out_)
    : in(in_)
    , out(out_) {}

Algorithm::~Algorithm() {}

/**
 * \brief The dryRun estimates resource consumption, especially
 * memory and processor time.
 */
void Algorithm::dryRun() {
  LOG(0, getName()) << "dry run not implemented" << std::endl;
}

AlgorithmFactory::AlgorithmMap *AlgorithmFactory::algorithmMap;
